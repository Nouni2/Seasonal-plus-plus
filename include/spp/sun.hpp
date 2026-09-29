#pragma once

#include "spp/core/time.hpp"
#include "spp/core/types.hpp"
#include "spp/earth.hpp"

// The Sun: its geocentric position, and its place in an observer's sky.
namespace spp::sun {

// Geometric geocentric ecliptic position (VSOP87D Earth, inverted), referred
// to the FK5 system.
Ecliptic geometric(Instant<TT> t);

// Apparent ecliptic position: geometric + nutation in longitude - aberration.
Ecliptic apparent_ecliptic(Instant<TT> t, const earth::Nutation& n);

// Apparent geocentric equatorial position (true equator and equinox of date).
Equatorial apparent(Instant<TT> t, const earth::Nutation& n);
Equatorial apparent(Instant<TT> t);   // computes nutation itself

// Apparent geocentric semidiameter of the Sun: 959.63" / R.
double semidiameter_deg(Instant<TT> t);

// The Sun in the observer's sky: topocentric azimuth and refracted altitude.
Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm = {});

}  // namespace spp::sun
