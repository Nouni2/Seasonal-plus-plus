#pragma once

#include "vector.hpp"
#include "spp/core/constants.hpp"

namespace spp::detail {
inline constexpr double light_speed_km_day=speed_of_light_m_s*seconds_per_day/1000;

// Relativistic aberration of a source direction for a moving receiver.
// ray must be nonzero and the receiver velocity (km/day) subluminal.
// The length of the ray is retained.
inline Vector aberration(Vector ray,Vector velocity_km_day) {
    const double range=norm(ray);
    const Vector n=(1.0/range)*ray;
    const Vector beta=(1.0/light_speed_km_day)*velocity_km_day;
    const double b=std::sqrt(1-dot(beta,beta)),nb=dot(n,beta);
    return (range/(1+nb))*(b*n+(1+nb/(1+b))*beta);
}

// Tangent differential on the unit sphere; n/apparent are the directions
// before/after the same aberration map, and dot(n,tangent)==0.
inline Vector aberration_tangent(Vector n,Vector apparent,Vector tangent,Vector velocity_km_day) {
    const Vector beta=(1.0/light_speed_km_day)*velocity_km_day;
    const double b=std::sqrt(1-dot(beta,beta)),u=dot(tangent,beta);
    return (1.0/(1+dot(n,beta)))*(b*tangent+(u/(1+b))*beta-u*apparent);
}

}  // namespace spp::detail
