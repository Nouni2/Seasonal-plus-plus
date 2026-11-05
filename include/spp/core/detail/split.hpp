#pragma once

#include <cmath>
#include <cstdint>

#include "spp/core/time.hpp"

namespace spp::detail {

inline constexpr double kMaxSplitDays = 0x1p62;

// Shared by the compiled pipeline and header-only ephemeris adapters. Leave
// unnormalizable splits unchanged so the caller can reject them without UB.
inline void normalize(std::int64_t& day, double& frac) {
    if (frac >= 0.0 && frac < 1.0) return;
    if (!(std::fabs(frac) < kMaxSplitDays) || day > (std::int64_t{1} << 62)
        || day < -(std::int64_t{1} << 62)) return;
    const double whole = std::floor(frac);
    day += static_cast<std::int64_t>(whole);
    frac -= whole;
    if (frac >= 1.0) { day += 1; frac -= 1.0; }
}

template<class Scale>
Instant<Scale> normalized(Instant<Scale> t) {
    normalize(t.day, t.frac);
    return t;
}

}  // namespace spp::detail
