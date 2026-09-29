#pragma once

#include <cmath>

#include "spp/core/constants.hpp"

namespace spp {

// Normalize to [0, 360).
inline double norm360(double deg) {
    double r = std::fmod(deg, 360.0);
    return r < 0.0 ? r + 360.0 : r;
}

// Wrap to [-180, 180).
inline double wrap180(double deg) {
    return norm360(deg + 180.0) - 180.0;
}

}  // namespace spp
