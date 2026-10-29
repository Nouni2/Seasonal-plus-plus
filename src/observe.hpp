#pragma once

// Internal: the position pipeline, in one place. sun::horizontal and the event
// solvers both observe the Sun through these two functions.

#include "spp/core/time.hpp"
#include "spp/core/types.hpp"
#include "spp/earth.hpp"

namespace spp::detail {

// Everything about one instant that does not depend on the observer. Nutation
// is computed once and shared by the apparent place and the equation of the
// equinoxes, so both always use the same Δψ.
struct Geocentric {
    earth::Nutation nutation;
    Equatorial sun;      // apparent geocentric place, true equator and equinox of date
    Instant<UT1> ut1;
};

Geocentric geocentric(const Time& t);

// Local apparent sidereal time at the instant of `g`, degrees.
double lst_deg(const Geocentric& g, double lon_east_deg);

// What the observer sees: topocentric azimuth and altitude, refracted unless
// atm.refraction is false.
Horizontal observed(const Geocentric& g, const Observer& obs, const Atmosphere& atm);

}  // namespace spp::detail
