#include "spp/model/espenak_meeus.hpp"

#include <array>
#include <cmath>
#include <cstddef>

namespace spp::model::espenak_meeus {

namespace {

// Observed Delta T (TT - UT1) at the start of each year, seconds (IERS / USNO).
// 2026: 69.184 s - (UT1 - UTC), with UT1 - UTC = +0.07 s from IERS Bulletin A
// (January 2026).
constexpr int kTableStartYear = 2005;
constexpr std::array kObserved = {
    64.69, 64.85, 65.15, 65.46, 65.78,   // 2005-2009
    66.07, 66.32, 66.60, 66.91, 67.28,   // 2010-2014
    67.64, 68.10, 68.59, 68.97, 69.22,   // 2015-2019
    69.36, 69.36, 69.29, 69.18, 69.10,   // 2020-2024
    69.07, 69.11,                        // 2025, 2026
};
constexpr int kTableEndYear = kTableStartYear + static_cast<int>(kObserved.size()) - 1;

double p(double x, double n) { return std::pow(x, n); }

}  // namespace

double delta_t_seconds(double y) {
    // Observed era: linear interpolation of the annual values (< ~0.2 s error).
    if (kTableStartYear <= y && y < kTableEndYear) {
        const auto idx = static_cast<std::size_t>(y - kTableStartYear);
        const double frac = (y - kTableStartYear) - static_cast<double>(idx);
        const double v0 = kObserved[idx];
        const double v1 = kObserved[idx + 1];
        return v0 + frac * (v1 - v0);
    }

    // Near future: from the last observation with zero slope (observed ΔT has
    // been flat since 2021), curving up to the Espenak-Meeus 93 s at 2050 so the
    // next piece stays continuous. A prediction: the uncertainty grows with time.
    if (kTableEndYear <= y && y < 2050) {
        const double last = kObserved.back();
        const double u = (y - kTableEndYear) / (2050.0 - kTableEndYear);
        return last + (93.0 - last) * u * u;
    }

    // Prediction era: parabola plus the continuity correction at 2050.
    if (2050 <= y && y < 2150) {
        const double t = (y - 1820.0) / 100.0;
        return -20 + 32 * t * t - 0.5628 * (2150.0 - y);
    }

    if (1986 <= y && y < 2005) {
        const double t = y - 2000.0;
        return 63.86 + 0.3345 * t - 0.060374 * t * t + 0.0017275 * p(t, 3) + 0.000651814 * p(t, 4)
               + 0.00002373599 * p(t, 5);
    }

    // Historical records (Espenak & Meeus fits)
    if (1961 <= y && y < 1986) {
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

    // Pre-telescopic (continuous at 1600, 500 and -500)
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
