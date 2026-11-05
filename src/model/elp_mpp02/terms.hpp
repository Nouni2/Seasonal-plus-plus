#pragma once

#include <array>

namespace spp::model::elp_mpp02::detail {

struct MainTerm {
    std::array<int, 4> argument;
    // A, B1..B6 in the upstream units (rad for angles, km for distance).
    // The published adjustment uses B1..B5; B6 is retained for provenance.
    std::array<double, 7> coefficient;
};

struct PerturbationTerm {
    // D, F, l, l', Mercury, Venus, Earth-Moon, Mars, Jupiter, Saturn,
    // Uranus, Neptune, zeta.
    std::array<int, 13> argument;
    double amplitude;
    double phase_rad;
};

}  // namespace spp::model::elp_mpp02::detail
