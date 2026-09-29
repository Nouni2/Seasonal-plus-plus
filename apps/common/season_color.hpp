#pragma once

// Season colour shared by spp_spiral and spp_figures: blue around January,
// orange around July.

#include <array>
#include <cmath>
#include <numbers>

namespace common {

// Red, green, blue in [0, 1] for day-of-year `doy` (1-based) in a year of `days` days.
inline std::array<double, 3> season_rgb(int doy, int days) {
    const double s = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * (doy - 1) / days);
    constexpr std::array<double, 3> blue{0.30, 0.60, 1.00};
    constexpr std::array<double, 3> orange{1.00, 0.62, 0.20};
    return {blue[0] + (orange[0] - blue[0]) * s, blue[1] + (orange[1] - blue[1]) * s,
            blue[2] + (orange[2] - blue[2]) * s};
}

}  // namespace common
