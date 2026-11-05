#pragma once

#include "../core/gravity.hpp"
#include "spp/core/constants.hpp"

namespace spp::model::de440_solar {
using spp::detail::Vector;
using namespace spp::detail;
inline constexpr double radius_km=696000,j2=2.1961391516529825e-7;
inline constexpr double inertia_fraction=0.06884,spin_rad_day=14.1844*deg2rad;
inline Vector pole() {
    constexpr double ra=286.13*deg2rad,dec=63.87*deg2rad;
    return {std::cos(dec)*std::cos(ra),std::cos(dec)*std::sin(ra),std::sin(dec)};
}
// r and v are Sun-to-body ICRF vectors in km and km/TDB day.
inline Vector oblateness(Vector r,double gm_km3_day2) {
    constexpr std::array<double,3> zonal{0,0,j2};
    return zonal_acceleration(r,pole(),radius_km,gm_km3_day2,2,zonal);
}
inline Vector radiation(Vector r,double gm_km3_day2,double fraction) {
    // Repulsive photon pressure. The coefficient is a fraction of the
    // central gravitational magnitude, not an added attractive mass.
    const double r2=dot(r,r);
    return (fraction*gm_km3_day2/(r2*std::sqrt(r2)))*r;
}
inline Vector frame_dragging(Vector r,Vector v,double gm_km3_day2) {
    constexpr double c=299792.458*seconds_per_day;
    const auto angular_momentum=(gm_km3_day2*radius_km*radius_km*inertia_fraction*spin_rad_day)*pole();
    const double distance=norm(r);
    const auto omega=(1/(c*c*distance*distance*distance))
        *((-1.0)*angular_momentum+(3*dot(angular_momentum,r)/(distance*distance))*r);
    return 2*cross(omega,v);
}
} // namespace spp::model::de440_solar
