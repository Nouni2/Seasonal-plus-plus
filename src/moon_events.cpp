#include "spp/moon.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

#include "core/root_search.hpp"
#include "moon_internal.hpp"
#include "spp/core/angles.hpp"
#include "spp/sky.hpp"

namespace spp::moon::detail {
namespace {
constexpr double curvature=128; // |d² sin(altitude)/d(day)²|, terrestrial Moon
// Coordinate roundoff includes subtracting large heliocentric Earth vectors
// in the light-time calculation. This is separate from the requested angular
// contact tolerance. Keep the pure scalar search's smaller default budget.
constexpr double value_roundoff=1e-11;

double elapsed(const Time& a,const Time& b) {
    const auto x=a.utc(),y=b.utc();return static_cast<double>(x.day-y.day)+(x.frac-y.frac);
}

double geometric_target(const Atmosphere& effective,double target) {
    if(!effective.refraction || *effective.pressure_mbar==0)return target;
    const double scale=(*effective.pressure_mbar/1010)*(283/(273+effective.temp_c));
    if(!std::isfinite(scale) || scale>2)
        throw std::invalid_argument("Moon horizon search requires a monotonic terrestrial refraction configuration (scale <= 2)");
    auto f=[&](double h){return sky::refraction(h,effective)-target;};
    return spp::detail::refine_root(f,90,f(90),-90,f(-90),1e-12,100);
}
}  // namespace

HorizonResult horizon(const Time& start,const Time& end,const Observer& obs,
                      const Atmosphere& atm,const HorizonOptions& options,ModelRef model) {
    validate(obs,atm);
    if(!std::isfinite(options.altitude_deg) || std::abs(options.altitude_deg)>80
        || !std::isfinite(options.time_tolerance_s) || options.time_tolerance_s<1e-5 || options.time_tolerance_s>60
        || !std::isfinite(options.contact_tolerance_arcsec) || options.contact_tolerance_arcsec<1e-5
        || options.contact_tolerance_arcsec>1)
        throw std::invalid_argument("Moon horizon search: invalid altitude, time tolerance or contact tolerance");
    const auto a=start.utc(),b=end.utc();
    if(b.day<a.day || (b.day==a.day && b.frac<a.frac))
        throw std::invalid_argument("Moon horizon search: reversed interval");
    HorizonResult result;
    if(a.day==b.day && a.frac==b.frac)return result;
    (void)model(start.tdb());(void)model(end.tdb());
    const auto effective=sky::at_elevation(atm,obs.elev_m);
    const double target=geometric_target(effective,options.altitude_deg);
    const double tol=options.time_tolerance_s/seconds_per_day;
    const double contact=options.contact_tolerance_arcsec/3600*deg2rad
        *std::cos((std::abs(target)+0.5)*deg2rad);
    const Atmosphere airless{.refraction=false};
    auto condition=[&](const Coordinates& c) {
        const double limb=options.upper_limb?c.semidiameter_deg:0;
        return std::sin(c.geometric_alt_deg*deg2rad)-std::sin((target-limb)*deg2rad);
    };
    Time current=start;
    bool first=true;
    while(true) {
        const auto u=current.utc();
        const double remaining=b.day>u.day+1?2:static_cast<double>(b.day-u.day)+(b.frac-u.frac);
        const double width=std::min(1.0,remaining);
        std::unordered_map<double,Coordinates> memo;
        memo.reserve(128);
        auto at=[&](double x)->const Coordinates& {
            if(const auto found=memo.find(x);found!=memo.end())return found->second;
            auto place=coordinates(current.shifted(x),obs,airless,model);
            // Curvature/upper-limb bounds are specific to the terrestrial Moon,
            // not arbitrary near-Earth satellites supplied through ModelRef.
            if(place.distance_km<250000 || place.distance_km>600000)
                throw std::runtime_error("Moon horizon search: provider outside terrestrial lunar geometry");
            return memo.emplace(x,place).first->second;
        };
        auto append=[&](double offset,HorizonKind kind) {
            if(offset<0 || offset>=width)return;
            const auto& c=at(offset);
            const auto time=current.shifted(offset);
            const auto instant=time.utc();
            if(instant.day>b.day || (instant.day==b.day && instant.frac>=b.frac))return;
            result.events.push_back({kind,time,sky::refraction(c.geometric_alt_deg,effective),c.horizontal.az_deg});
        };
        const auto crossings=spp::detail::bounded_roots([&](double x){return condition(at(x));},width,tol,contact,curvature,value_roundoff);
        if(first) {
            // A tangent at start does not switch visibility. In particular,
            // a minimum exactly on the target is above immediately after
            // start, even though f(start)==0. Prefer the first crossing's
            // direction; with no crossing, sample away from any contact.
            const auto cross=std::find_if(crossings.begin(),crossings.end(),[](const auto& r){return !r.contact;});
            if(cross!=crossings.end())result.above_at_start=cross->x==0?cross->rising:!cross->rising;
            else {
                double value=condition(at(0));
                for(double x:{width/2,width}) {
                    const double candidate=condition(at(x));
                    if(std::abs(candidate)>std::abs(value))value=candidate;
                }
                result.above_at_start=value>0;
            }
            first=false;
        }
        for(const auto& root:crossings)append(root.x,root.contact?HorizonKind::Graze:
            root.rising?HorizonKind::Rise:HorizonKind::Set);
        const auto meridians=spp::detail::bounded_roots([&](double x){return std::sin(at(x).hour_angle_deg*deg2rad);},
            width,tol,1e-12,curvature,value_roundoff);
        for(const auto& root:meridians) {
            if(root.contact)throw std::runtime_error("Moon horizon search: stationary meridian alignment violates lunar motion assumptions");
            append(root.x,std::cos(at(root.x).hour_angle_deg*deg2rad)>0?HorizonKind::UpperTransit:HorizonKind::LowerTransit);
        }
        if(remaining<=1)break;
        current=current.shifted(1);
    }
    std::sort(result.events.begin(),result.events.end(),[](const auto& x,const auto& y){return elapsed(x.time,y.time)<0;});
    std::vector<HorizonEvent> unique;
    for(const auto& e:result.events) {
        bool duplicate=false;
        for(auto i=unique.rbegin();i!=unique.rend() && elapsed(e.time,i->time)<=tol;++i)
            if(i->kind==e.kind){duplicate=true;break;}
        if(!duplicate)unique.push_back(e);
    }
    result.events=std::move(unique);
    bool up=result.above_at_start;
    Time cursor=start;
    for(const auto& e:result.events)if(e.kind==HorizonKind::Rise || e.kind==HorizonKind::Set) {
        if(up)result.duration_hours+=elapsed(e.time,cursor)*24;
        up=e.kind==HorizonKind::Rise;
        if(elapsed(e.time,start)==0)result.above_at_start=up;
        cursor=e.time;
    }
    if(up)result.duration_hours+=elapsed(end,cursor)*24;
    return result;
}
}  // namespace spp::moon::detail
