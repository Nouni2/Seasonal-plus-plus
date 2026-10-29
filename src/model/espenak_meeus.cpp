#include "spp/model/espenak_meeus.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "delta_t_observed.hpp"

namespace spp::model::espenak_meeus {

namespace {

using detail::kObservedDeltaT;

constexpr int kTableFirstMjd = kObservedDeltaT[0].mjd;
constexpr int kTableLastMjd = std::end(kObservedDeltaT)[-1].mjd;

// Decimal year y = 2000 + (JD - 2451545) / 365.25 and MJD = JD - 2400000.5.
constexpr double kMjdJ2000 = 51544.5;
constexpr double mjd_of(double y) { return kMjdJ2000 + (y - 2000.0) * 365.25; }
constexpr double year_of(double mjd) { return 2000.0 + (mjd - kMjdJ2000) / 365.25; }

constexpr double kTableFirstYear = year_of(kTableFirstMjd);
constexpr double kTableLastYear = year_of(kTableLastMjd);

// Linear interpolation of the monthly table. The year -> MJD conversion can
// round just outside the table at its ends; the end intervals then extrapolate
// by that rounding amount.
double observed(double mjd) {
    const auto* first = std::begin(kObservedDeltaT);
    const auto* last = std::end(kObservedDeltaT) - 1;
    const auto* hi = std::upper_bound(first, last + 1, mjd,
                                      [](double m, const detail::ObservedDeltaT& e) { return m < e.mjd; });
    hi = std::clamp(hi, first + 1, last);
    const auto* lo = hi - 1;
    const double frac = (mjd - lo->mjd) / static_cast<double>(hi->mjd - lo->mjd);
    return lo->seconds + frac * (hi->seconds - lo->seconds);
}

double p(double x, double n) { return std::pow(x, n); }

}  // namespace

double delta_t_seconds(double y) {
    // Observed era: monthly IERS values, interpolated linearly. The month
    // after the last entry is not covered: the bridge below starts there.
    if (kTableFirstYear <= y && y < kTableLastYear) return observed(mjd_of(y));

    // Near future: from the last observation with zero slope, curving up to
    // the Espenak-Meeus 93 s at 2050 so the next piece stays continuous. A
    // prediction: the uncertainty grows with time.
    if (kTableLastYear <= y && y < 2050) {
        const double last = std::end(kObservedDeltaT)[-1].seconds;
        const double u = (y - kTableLastYear) / (2050.0 - kTableLastYear);
        return last + (93.0 - last) * u * u;
    }

    // Prediction era: parabola plus the continuity correction at 2050.
    if (2050 <= y && y < 2150) {
        const double t = (y - 1820.0) / 100.0;
        return -20 + 32 * t * t - 0.5628 * (2150.0 - y);
    }

    // Historical records (Espenak & Meeus fits)
    if (1961 <= y && y < kTableFirstYear) {
        const double t = y - 1975.0;
        return 45.45 + 1.067 * t - p(t, 2) / 260.0 - p(t, 3) / 718.0;
    }
    if (1941 <= y && y < 1961) {
        const double t = y - 1950.0;
        return 29.07 + 0.407 * t - p(t, 2) / 233.0 + p(t, 3) / 2547.0;
    }
    if (1920 <= y && y < 1941) {
        const double t = y - 1920.0;
        return 21.20 + 0.84493 * t - 0.076100 * p(t, 2) + 0.0020936 * p(t, 3);
    }
    if (1900 <= y && y < 1920) {
        const double t = y - 1900.0;
        return -2.79 + 1.494119 * t - 0.0598939 * p(t, 2) + 0.0061966 * p(t, 3) - 0.000197 * p(t, 4);
    }
    if (1860 <= y && y < 1900) {
        const double t = y - 1860.0;
        return 7.62 + 0.5737 * t - 0.251754 * p(t, 2) + 0.01680668 * p(t, 3) - 0.0004473624 * p(t, 4)
               + p(t, 5) / 233174.0;
    }
    if (1800 <= y && y < 1860) {
        // The t^7 term is essential: without it the fit diverges by 1860.
        const double t = y - 1800.0;
        return 13.72 - 0.332447 * t + 0.0068612 * p(t, 2) + 0.0041116 * p(t, 3) - 0.00037436 * p(t, 4)
               + 0.0000121272 * p(t, 5) - 0.0000001699 * p(t, 6) + 0.000000000875 * p(t, 7);
    }
    if (1700 <= y && y < 1800) {
        const double t = y - 1700.0;
        return 8.83 + 0.1603 * t - 0.0059285 * p(t, 2) + 0.00013336 * p(t, 3) - p(t, 4) / 1174000.0;
    }
    if (1600 <= y && y < 1700) {
        const double t = y - 1600.0;
        return 120.0 - 0.9808 * t - 0.01532 * p(t, 2) + p(t, 3) / 7129.0;
    }

    // Pre-telescopic. The published pieces do not join exactly: ΔT jumps by up
    // to 0.25 s at the piece boundaries (README §2.6, tools/reference_check.py).
    if (500 <= y && y < 1600) {
        const double u = (y - 1000.0) / 100.0;
        return 1574.2 - 556.01 * u + 71.23472 * p(u, 2) + 0.319781 * p(u, 3) - 0.8503463 * p(u, 4)
               - 0.005050998 * p(u, 5) + 0.0083572073 * p(u, 6);
    }
    if (-500 <= y && y < 500) {
        const double u = y / 100.0;
        return 10583.6 - 1014.41 * u + 33.78311 * p(u, 2) - 5.952053 * p(u, 3) - 0.1798452 * p(u, 4)
               + 0.022174192 * p(u, 5) + 0.0090316521 * p(u, 6);
    }

    // Before -500 and from 2150 on: long-term tidal parabola.
    if (y < -500 || y >= 2150) {
        const double t = (y - 1820.0) / 100.0;
        return -20 + 32 * t * t;
    }
    return 0.0;   // only reached for NaN
}

}  // namespace spp::model::espenak_meeus
