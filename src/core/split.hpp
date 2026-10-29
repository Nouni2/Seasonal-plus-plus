#pragma once

// Internal: normalization of split Julian Dates.

#include <cmath>
#include <cstdint>

#include "spp/core/time.hpp"

namespace spp::detail {

// Days beyond which a split cannot be normalized without integer overflow.
inline constexpr double kMaxSplitDays = 0x1p62;

// Brings frac into [0, 1), carrying whole days into day. frac - floor(frac) is
// exact; the one rounding case is a tiny negative frac, whose sum with 1.0
// rounds to 1.0 and is carried as a whole day. A non-finite frac, or one or a
// day too large to carry without overflow, is left as it is: the results are
// then NaN or meaningless, but defined.
inline void normalize(std::int64_t& day, double& frac) {
    if (frac >= 0.0 && frac < 1.0) return;   // already normalized: the common case
    if (!(std::fabs(frac) < kMaxSplitDays) || day > (std::int64_t{1} << 62) || day < -(std::int64_t{1} << 62))
        return;
    const double whole = std::floor(frac);
    day += static_cast<std::int64_t>(whole);
    frac -= whole;
    if (frac >= 1.0) {
        day += 1;
        frac -= 1.0;
    }
}

// Instant is an aggregate, so it can hold any frac. Every function that takes an
// Instant normalizes it first, so all splits of one instant give the same result.
// A normalized Instant is returned unchanged.
template <class Scale>
Instant<Scale> normalized(Instant<Scale> t) {
    normalize(t.day, t.frac);
    return t;
}

}  // namespace spp::detail
