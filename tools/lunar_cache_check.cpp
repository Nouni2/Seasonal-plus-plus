// Interpolation error is checked separately from the ELP/JPL model residuals.
#include "spp/spp.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "error_stats.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <format>
#include <limits>
#include <random>
#include <stdexcept>
#include <type_traits>

namespace {
using spp::model::elp_mpp02::CachedEphemeris;
using spp::model::elp_mpp02::Ephemeris;
using spp::model::elp_mpp02::Fit;
using Instant = spp::Instant<spp::TDB>;
static_assert(spp::moon::Ephemeris<CachedEphemeris>);
static_assert(spp::moon::PositionEphemeris<CachedEphemeris>);
static_assert(!std::is_invocable_v<CachedEphemeris&, spp::Instant<spp::TT>>);

struct OtherFrame {
    spp::CartesianState<void,void> operator()(Instant);
};
struct OtherScale {
    spp::moon::State operator()(spp::Instant<spp::TT>);
};
using OtherScaleCache=spp::model::HermiteEphemeris<OtherScale,spp::TT>;
static_assert(!spp::moon::Ephemeris<spp::model::HermiteEphemeris<OtherFrame>>);
static_assert(std::is_invocable_v<OtherScaleCache&,spp::Instant<spp::TT>>);
static_assert(!spp::moon::PositionSource<OtherScaleCache>);

struct NonFinite {
    spp::moon::State operator()(Instant) {
        return {{1,2,3},{0,std::numeric_limits<double>::infinity(),0}};
    }
};

void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}
template<class Exception, class F> void throws(F f) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error("expected cache exception not raised");
}
Instant shifted(Instant t, double days) {
    return spp::detail::normalized(Instant{t.day, t.frac+days});
}
double seconds(const spp::Time& a, const spp::Time& b) {
    const auto x=a.utc(), y=b.utc();
    return (static_cast<double>(x.day-y.day)+(x.frac-y.frac))*spp::seconds_per_day;
}

// An exact cubic checks interpolation and its derivative without sharing the
// production implementation's coefficient construction or lunar equations.
struct Polynomial {
    int* calls;
    bool* fail;
    spp::moon::State operator()(Instant t) {
        ++*calls;
        if (t.day<2451545 || (t.day==2451545 && t.frac<0.12517)
            || t.day>2451546 || (t.day==2451546 && t.frac>0.37543))
            throw std::out_of_range("synthetic provider sampled outside its domain");
        if (*fail && *calls%2==0) throw std::runtime_error("injected endpoint failure");
        const double u=static_cast<double>(t.day-2451545)+t.frac;
        return {{300000+100*u+20*u*u+u*u*u,
                 -200000+200*u-10*u*u+2*u*u*u,
                 1000+30*u+5*u*u-3*u*u*u},
                {100+40*u+3*u*u, 200-20*u+6*u*u, 30+10*u-9*u*u}};
    }
};

void contracts() {
    int calls=0, expected_calls=0;
    bool fail=false;
    const Instant first{2451545,0.12517}, last{2451546,0.37543};
    spp::model::HermiteEphemeris<Polynomial> cached(Polynomial{&calls,&fail},first,last);
    Polynomial direct{&expected_calls,&fail};
    for (const Instant t : {first, shifted(first,1e-8), Instant{2451545,0.25},
                            Instant{2451545,std::nextafter(1.0,0.0)},
                            Instant{2451546,0.0}, shifted(last,-1e-8), last}) {
        const auto a=cached(t), b=direct(t);
        for (std::size_t q=0;q<3;++q) {
            require(std::abs(a.position_km[q]-b.position_km[q])<1e-8,"cubic position/endpoint mismatch");
            require(std::abs(a.velocity_km_per_day[q]-b.velocity_km_per_day[q])<1e-6,
                    "cubic derivative/endpoint mismatch");
        }
    }
    cached.clear(); calls=0;
    (void)cached.position({2451545,0.2501});
    require(calls==2,"cache miss did not use exactly two state samples");
    (void)cached.position({2451545,0.2502});
    (void)cached({2451545,0.2503});
    (void)cached.position({2451544,1.2502});
    require(calls==2,"cache hit or normalized split resampled the provider");
    (void)cached.position({2451545,0.254});
    require(calls==3,"adjacent interval did not reuse its common endpoint");
    cached.clear();
    (void)cached.position({2451545,0.2502});
    require(calls==5,"clear retained a populated segment or endpoint");
    cached.clear();
    for(int i=0;i<256;++i)(void)cached.position({2451545,0.1252+i*0.0035});
    const int populated=calls;
    for(int i=0;i<256;++i)(void)cached.position({2451545,0.1252+i*0.0035});
    require(calls==populated,"hash collisions evicted intervals that fit the cache capacity");
    for(int round=0;round<5;++round)for(int i=0;i<1000;++i)
        (void)cached.position({2451545,0.1252+static_cast<double>(i)/1000});
    const int before_eviction=calls;
    const auto evicted=cached({2451545,0.1252});
    const auto exact=direct({2451545,0.1252});
    require(calls>before_eviction,"eviction/reindex check did not refill an expired interval");
    require(std::abs(evicted.position_km[0]-exact.position_km[0])<1e-8,
            "cache eviction returned another segment's state");
    cached.clear(); calls=0; fail=true;
    throws<std::runtime_error>([&]{cached.position({2451545,0.2502});});
    fail=false;
    (void)cached.position({2451545,0.2502});
    require(calls==4,"endpoint failure published a partially built segment");

    for (double frac : {std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()})
        throws<std::invalid_argument>([&]{cached.position({2451545,frac});});
    for (Instant t : {shifted(first,-1e-8),shifted(last,1e-8),
                       Instant{std::numeric_limits<std::int64_t>::max(),0},
                       Instant{std::numeric_limits<std::int64_t>::min(),0},Instant{2451545,1e100}})
        throws<std::out_of_range>([&]{cached.position(t);});
    spp::model::HermiteEphemeris<NonFinite> invalid(NonFinite{},first,last);
    throws<std::invalid_argument>([&]{invalid.position(first);});
    throws<std::invalid_argument>([]{CachedEphemeris invalid_fit(static_cast<Fit>(42));});

    CachedEphemeris de405, llr(Fit::LLR);
    const Instant epoch{2451545,0.5001};
    require(de405.position(epoch).position_km!=llr.position(epoch).position_km,"cache mixed model fits");
    (void)de405.position(shifted(epoch,-1.5/spp::seconds_per_day));
    (void)de405.position(shifted(epoch,-1.4/spp::seconds_per_day));
    require(de405.source_evaluations()==2,"light-time calls rebuilt an already sampled interval");
    for (Instant t : {Ephemeris::first_epoch,shifted(Ephemeris::first_epoch,1e-10),
                       shifted(Ephemeris::last_epoch,-1e-10),Ephemeris::last_epoch})
        (void)de405(t);
    throws<std::out_of_range>([&]{de405.position(shifted(Ephemeris::first_epoch,-1e-8));});
    throws<std::out_of_range>([&]{de405.position(shifted(Ephemeris::last_epoch,1e-8));});
    std::puts("Cache contracts: cubic/derivative, split dates, bounds, reuse, eviction, fits, exceptions: PASS");
}

void interpolation() {
    spp::check::ErrorStats modern_p, modern_v, historical_p, historical_v;
    std::mt19937_64 rng(49);
    for (Fit fit : {Fit::DE405,Fit::LLR}) {
        const Ephemeris direct(fit);
        CachedEphemeris cached(fit);
        auto check=[&](Instant t) {
            const auto a=cached(t), b=direct(t);
            const auto p=cached.position(t);
            const auto label=std::format("{} TDB {}+{:.15f}",fit==Fit::LLR?"LLR":"DE405",t.day,t.frac);
            const bool modern=t.day>=2451545 && t.day<=2488070;
            for (std::size_t q=0;q<3;++q) {
                require(a.position_km[q]==p.position_km[q],"cached position/state paths differ");
                const auto component=std::format("{} {}",label,"xyz"[q]);
                (modern?modern_p:historical_p).add((a.position_km[q]-b.position_km[q])*1000,component);
                (modern?modern_v:historical_v).add(a.velocity_km_per_day[q]-b.velocity_km_per_day[q],component);
            }
        };
        for (int i=0;i<192;++i) {
            const bool modern=i<64;
            const auto first=modern?2451545LL:Ephemeris::first_epoch.day;
            const auto span=modern?36524LL:Ephemeris::last_epoch.day-first;
            const auto day=first+static_cast<std::int64_t>(rng()%static_cast<std::uint64_t>(span));
            constexpr int intervals=CachedEphemeris::intervals_per_day;
            const int bin=static_cast<int>(rng()%intervals);
            // Interior points, endpoints, both sides of a segment boundary.
            for (double u : {0.0,0.125,0.5,0.875,std::nextafter(1.0,0.0)})
                check(spp::detail::normalized(Instant{day,(bin+u)/intervals}));
            const double edge=static_cast<double>(bin+1)/intervals;
            check({day,std::nextafter(edge,0.0)});
            check(spp::detail::normalized(Instant{day,edge}));
            if (edge<1.0 || day+1<Ephemeris::last_epoch.day)
                check(spp::detail::normalized(Instant{day,std::nextafter(edge,2.0)}));
        }
        // Sweep through a lunation with different offsets near daily boundaries.
        for (int i=0;i<128;++i)
            check({2461347+i/4,(i%4)*0.25+static_cast<double>(i%31)/86400.0});
        // Retain the difficult positions/rates exposed by the first five-minute
        // cache run when changing the production interpolation interval.
        for(Instant t : {Instant{2485482,0.782986111111111},Instant{2488048,0.336371527777778},
                         Instant{1402526,0.550347222222222},Instant{1518288,0.824652777777778}})check(t);
        for (Instant t : {Ephemeris::first_epoch,shifted(Ephemeris::first_epoch,1e-9),
                           shifted(Ephemeris::last_epoch,-1e-9),Ephemeris::last_epoch}) check(t);
    }
    std::puts("Interpolation residuals against direct full ELP (both fits; per Cartesian component):");
    bool pass=modern_p.report("2000-2100 position","m",0.001);
    pass=modern_v.report("2000-2100 velocity","km/day",1e-4)&&pass;
    // At distant epochs, rounding of ELP's century/argument polynomials also
    // affects the direct comparator. Keep that separate from the modern set.
    pass=historical_p.report("full-domain position","m",0.1)&&pass;
    pass=historical_v.report("full-domain velocity","km/day",0.02)&&pass;
    require(pass,"lunar interpolation exceeded its regression error budget");
}

void observations() {
    spp::check::ErrorStats sky, range, cycle, illumination, events, disk_angles, disk_light, limb;
    for (Fit fit : {Fit::DE405,Fit::LLR}) {
        const Ephemeris direct(fit);
        CachedEphemeris cached(fit);
        for (int i=0;i<12;++i) {
            const auto t=spp::Time::from_utc(2000+i*8,1+i%12,1+i%27,12,0,0.0,0.0);
            for (spp::Observer obs : {spp::Observer{48.8566,2.3522,35},spp::Observer{-78,-170,5000}}) {
                const auto a=spp::moon::observe(t,obs,cached,{.refraction=false});
                const auto b=spp::moon::observe(t,obs,direct,{.refraction=false});
                const auto h=spp::moon::horizontal(t,obs,cached,{.refraction=false});
                require(h.alt_deg==a.horizontal.alt_deg && h.az_deg==a.horizontal.az_deg,
                        "cache changed observer-path consistency");
                const auto label=std::format("{} {} lat {}",fit==Fit::LLR?"LLR":"DE405",t.iso_utc(),obs.lat_deg);
                sky.add(3600*std::hypot(a.horizontal.alt_deg-b.horizontal.alt_deg,
                    spp::wrap180(a.horizontal.az_deg-b.horizontal.az_deg)*std::cos(b.horizontal.alt_deg*spp::deg2rad)),label);
                range.add((a.distance_km-b.distance_km)*1000,label);
                cycle.add(spp::wrap180(a.phase.cycle_deg-b.phase.cycle_deg)*3600,label);
                illumination.add((a.phase.illuminated_fraction-b.phase.illuminated_fraction)*100,label);
                disk_angles.add((a.disk.phase_angle_deg-b.disk.phase_angle_deg)*3600,label);
                disk_light.add((a.disk.illuminated_fraction-b.disk.illuminated_fraction)*100,label);
                if(a.disk.bright_limb_defined && b.disk.bright_limb_defined)
                    limb.add(spp::wrap180(a.disk.bright_limb_deg-b.disk.bright_limb_deg)*3600,label);
                // Independent Cartesian tangent-plane construction verifies
                // orientation/sign and the phase triangle used by the shader.
                const auto tt=t.tt();const auto nut=spp::earth::nutation(tt);
                const double lst=spp::earth::lst_deg(t.ut1(),nut,spp::norm360(obs.lon_deg));
                const auto solar=spp::sky::diurnal_aberration(spp::sky::parallax(spp::sun::apparent(tt,nut),obs,lst),obs,lst);
                auto unit=[](const spp::Equatorial& eq) {
                    const double r=eq.ra_deg*spp::deg2rad,d=eq.dec_deg*spp::deg2rad;
                    return std::array<double,3>{std::cos(d)*std::cos(r),std::cos(d)*std::sin(r),std::sin(d)};
                };
                auto dot=[](const auto& x,const auto& z){return x[0]*z[0]+x[1]*z[1]+x[2]*z[2];};
                const auto m=unit(a.equatorial),s=unit(solar);
                const double r=a.equatorial.ra_deg*spp::deg2rad,d=a.equatorial.dec_deg*spp::deg2rad;
                const std::array<double,3> east{-std::sin(r),std::cos(r),0},north{-std::sin(d)*std::cos(r),-std::sin(d)*std::sin(r),std::cos(d)};
                const double position_angle=std::atan2(dot(s,east),dot(s,north))*spp::rad2deg;
                require(std::abs(spp::wrap180(position_angle-a.disk.bright_limb_deg))<1e-9,"bright-limb projection/sign mismatch");
                std::array<double,3> sun_from_moon{};
                for(std::size_t q=0;q<3;++q)sun_from_moon[q]=solar.dist_au*spp::astronomical_unit_km*s[q]-a.distance_km*m[q];
                const double cos_phase=-dot(sun_from_moon,m)/std::sqrt(dot(sun_from_moon,sun_from_moon));
                // The apparent-direction triangle is only an approximation
                // to physical illumination: observer aberration must not
                // rotate the illumination geometry. Independent Horizons
                // Illu%/phi fixtures check the tighter physical result.
                require(std::abs(a.disk.illuminated_fraction-0.5*(1+cos_phase))<1e-4,"physical illumination far from its apparent-triangle approximation");
                const double phi=obs.lat_deg*spp::deg2rad,l=lst*spp::deg2rad;
                const std::array<double,3> zenith{std::cos(phi)*std::cos(l),std::cos(phi)*std::sin(l),std::sin(phi)};
                const double q=std::atan2(dot(zenith,east),dot(zenith,north))*spp::rad2deg;
                require(std::abs(spp::wrap180(q-a.disk.parallactic_angle_deg))<1e-9,"parallactic angle projection/sign mismatch");
            }
        }
        const auto start=spp::Time::from_utc(2026,11,1);
        const auto a=spp::moon::phases(start,start.shifted(35),cached);
        const auto b=spp::moon::phases(start,start.shifted(35),direct);
        require(a.size()==b.size() && !a.empty(),"cache changed principal-phase count");
        for (std::size_t i=0;i<a.size();++i) {
            require(a[i].phase==b[i].phase,"cache changed principal-phase ordering");
            events.add(seconds(a[i].time,b[i].time),a[i].time.iso_utc());
        }
    }
    std::puts("Added observation error from interpolation, separate from JPL/USNO model residuals:");
    bool pass=sky.report("sky direction","arcsec",1e-5);
    pass=range.report("observer range","m",0.001)&&pass;
    pass=cycle.report("phase cycle","arcsec",1e-5)&&pass;
    pass=illumination.report("illumination","percentage points",1e-7)&&pass;
    pass=disk_angles.report("disk phase angle","arcsec",1e-5)&&pass;
    pass=disk_light.report("disk illumination","percentage points",1e-7)&&pass;
    pass=limb.report("bright-limb direction","arcsec",1e-5)&&pass;
    pass=events.report("phase time vs direct","s",0.02)&&pass;
    require(pass,"interpolation changed lunar observations beyond its budget");
}
}  // namespace

int main() {
    bool failed=false;
    auto run=[&](const char* name,auto check) {
        try { check(); }
        catch (const std::exception& e) { std::fprintf(stderr,"FAIL [%s]: %s\n",name,e.what()); failed=true; }
    };
    run("cache contracts",contracts);
    run("interpolation",interpolation);
    run("cached observations",observations);
    if (failed) return 1;
    std::puts("PASS: interpolation budgets and cache contracts (sampled, not global error guarantees)");
    return 0;
}
