#pragma once

#include "spp/core/time.hpp"
#include "spp/core/types.hpp"
#include "spp/earth.hpp"

// The Sun as seen from the geocenter.
namespace spp::sun {

// Geometric geocentric ecliptic position (VSOP87D Earth, inverted).
Ecliptic geometric(Instant<TT> t);

// Apparent ecliptic position: geometric + nutation in longitude - aberration.
Ecliptic apparent_ecliptic(Instant<TT> t, const earth::Nutation& n);

// Apparent geocentric equatorial position (true equator and equinox of date).
Equatorial apparent(Instant<TT> t, const earth::Nutation& n);
Equatorial apparent(Instant<TT> t);   // computes nutation itself

}  // namespace spp::sun
