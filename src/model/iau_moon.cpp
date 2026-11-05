#include "spp/model/iau_moon.hpp"

#include <cmath>
#include <stdexcept>
#include "spp/core/angles.hpp"
#include "spp/core/detail/split.hpp"

namespace spp::model::iau_moon {
namespace {
// Published numerical data: NASA/NAIF pck00011.tpc, BODY301 and BODY3.
// Fundamental angular rates below are the original degrees/TDB day.
constexpr std::array<double,13> phase{125.045,250.089,260.008,176.625,357.529,
    311.589,134.963,276.617,34.226,15.134,119.743,239.961,25.053};
constexpr std::array<double,13> rate{-0.0529921,-0.1059842,13.0120009,13.3407154,0.9856003,
    26.4057084,13.0649930,0.3287146,1.7484877,-0.1589763,0.0036096,0.1643573,12.9590088};
constexpr std::array<double,13> ra{-3.8787,-0.1204,0.0700,-0.0172,0,0.0072,0,0,0,-0.0052,0,0,0.0043};
constexpr std::array<double,13> dec{1.5419,0.0239,-0.0278,0.0068,0,-0.0029,0.0009,0,0,0.0008,0,0,-0.0009};
constexpr std::array<double,13> pm{3.5610,0.1208,-0.0642,0.0158,0.0252,-0.0066,-0.0047,-0.0046,
    0.0028,0.0052,0.0040,0.0019,-0.0044};

double linear(double offset,double slope,std::int64_t days,double fraction) {
    // The float head times the bounded integer day is exact. Reduce that
    // product before adding the tail and fractional day to avoid losing
    // angular resolution in a many-million-degree prime meridian.
    const double head=static_cast<double>(static_cast<float>(slope)),tail=slope-head;
    return norm360(offset+std::fmod(head*static_cast<double>(days),360)
        +tail*static_cast<double>(days)+slope*fraction);
}
struct A {};struct B {}; // local intermediate rotation axes
template<class From,class To>
FrameRotation<From,To> z(double angle,double velocity) {
    const double c=std::cos(angle),s=std::sin(angle);
    return {{{{c,s,0},{-s,c,0},{0,0,1}}},
        {{{-s*velocity,c*velocity,0},{-c*velocity,-s*velocity,0},{0,0,0}}}};
}
template<class From,class To>
FrameRotation<From,To> x(double angle,double velocity) {
    const double c=std::cos(angle),s=std::sin(angle);
    return {{{{1,0,0},{0,c,s},{0,-s,c}}},
        {{{0,0,0},{0,-s*velocity,c*velocity},{0,-c*velocity,-s*velocity}}}};
}
}

Orientation Rotation::operator()(Instant<TDB> t) const {
    if(!std::isfinite(t.frac))throw std::invalid_argument("IAU Moon rotation: TDB fraction must be finite");
    t=spp::detail::normalized(t);
    if(t.frac<0 || t.frac>=1 || t.day<first_epoch.day || t.day>last_epoch.day
        || (t.day==last_epoch.day && t.frac>0))
        throw std::out_of_range("IAU Moon rotation: epoch outside -31..11 Julian centuries from J2000");
    const auto whole=t.day-jd_j2000_day;
    const double d=static_cast<double>(whole)+t.frac,T=d/days_per_century;
    Orientation out{269.9949+0.0031*T,66.5392+0.0130*T,
        linear(38.3213,13.17635815,whole,t.frac)-1.4e-12*d*d,
        0.0031/days_per_century,0.0130/days_per_century,13.17635815-2.8e-12*d,{}};
    for(std::size_t i=0;i<phase.size();++i) {
        const double angle=linear(phase[i],rate[i],whole,t.frac)*deg2rad;
        const double s=std::sin(angle),c=std::cos(angle),v=rate[i]*deg2rad;
        out.pole_ra_deg+=ra[i]*s;out.pole_dec_deg+=dec[i]*c;out.prime_meridian_deg+=pm[i]*s;
        out.pole_ra_rate_deg_day+=ra[i]*c*v;
        out.pole_dec_rate_deg_day-=dec[i]*s*v;
        out.prime_meridian_rate_deg_day+=pm[i]*c*v;
    }
    out.pole_ra_deg=norm360(out.pole_ra_deg);out.prime_meridian_deg=norm360(out.prime_meridian_deg);
    out.rotation=z<B,frame::MoonMeanEarth>(out.prime_meridian_deg*deg2rad,out.prime_meridian_rate_deg_day*deg2rad)
        *x<A,B>((90-out.pole_dec_deg)*deg2rad,-out.pole_dec_rate_deg_day*deg2rad)
        *z<frame::Icrf,A>((90+out.pole_ra_deg)*deg2rad,out.pole_ra_rate_deg_day*deg2rad);
    return out;
}
}  // namespace spp::model::iau_moon
