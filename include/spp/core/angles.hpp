#pragma once

#include <cmath>

#include "spp/core/constants.hpp"

namespace spp {

// Normalize to [0, 360). A negative remainder within half an ulp of 0 rounds
// to 360 when 360 is added, and is returned as 0; -0.0 is returned as +0.0;
// NaN stays NaN.
inline double norm360(double deg) {
    double r = std::fmod(deg, 360.0);
    if (r < 0.0) r += 360.0;
    return r >= 360.0 ? 0.0 : r + 0.0;
}

// Wrap to [-180, 180).
inline double wrap180(double deg) {
    return norm360(deg + 180.0) - 180.0;
}

}  // namespace spp
