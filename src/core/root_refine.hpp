#pragma once

#include <cmath>

namespace spp::detail {

// Illinois refinement. Endpoints are ordered by sign, not by time:
// f(lo) >= 0 > f(hi). Returns a zero or the final bracket midpoint.
template<class F>
double refine_root(F f, double lo, double f_lo, double hi, double f_hi, double tol_days, int max_iter) {
    double x = 0.5 * (lo + hi);
    int side = 0;
    for (int iter = 0; iter < max_iter && std::abs(hi - lo) > tol_days; ++iter) {
        x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if (!((lo < x && x < hi) || (hi < x && x < lo))) x = 0.5 * (lo + hi);
        const double fx = f(x);
        if (fx == 0.0) return x;
        if (fx > 0.0) {
            lo = x;
            f_lo = fx;
            if (side == 1) f_hi *= 0.5;
            side = 1;
        } else {
            hi = x;
            f_hi = fx;
            if (side == -1) f_lo *= 0.5;
            side = -1;
        }
    }
    return 0.5 * (lo + hi);
}

}  // namespace spp::detail
