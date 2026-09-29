#pragma once

// Internal: the expensive, observer-independent series of one instant.
//
// Direct evaluation costs ~2,500 cosines (VSOP87D) plus the nutation series.
// Both are smooth over days, so they are served from Chebyshev fits over
// short time segments, built on demand from the direct series and cached per
// thread. The fit error is kept at the level of double rounding in the direct
// sums (tools/cache_check.cpp measures it).

#include "spp/core/time.hpp"

namespace spp::detail {

struct EarthState {
    double l_rad;         // VSOP87D Earth longitude (not normalized)
    double b_rad;         // latitude
    double r_au;          // radius vector
    double dpsi_arcsec;   // IAU 1980 nutation in longitude
    double deps_arcsec;   // nutation in obliquity
};

// Direct series evaluation. T: Julian centuries of TT.
EarthState earth_state_direct(double T);

// Cached Chebyshev evaluation.
EarthState earth_state(Instant<TT> t);

}  // namespace spp::detail
