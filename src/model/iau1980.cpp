#include "spp/model/iau1980.hpp"

#include <cmath>

#include "spp/core/angles.hpp"
#include "spp/core/constants.hpp"
#include "iau1980_terms.hpp"

namespace spp::model::iau1980 {

Nutation nutation(double T) {
    const double T2 = T * T;
    const double T3 = T2 * T;

    // Fundamental luni-solar arguments (Meeus 22.1), degrees -> radians.
    const double D  = norm360(297.85036 + 445267.111480 * T - 0.0019142 * T2 + T3 / 189474.0) * deg2rad;
    const double M  = norm360(357.52772 + 35999.050340 * T - 0.0001603 * T2 - T3 / 300000.0) * deg2rad;
    const double Mp = norm360(134.96298 + 477198.867398 * T + 0.0086972 * T2 + T3 / 56250.0) * deg2rad;
    const double F  = norm360(93.27191 + 483202.017538 * T - 0.0036825 * T2 + T3 / 327270.0) * deg2rad;
    const double Om = norm360(125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000.0) * deg2rad;

    double d_psi = 0.0;
    double d_eps = 0.0;
    for (const Row& r : kTerms) {
        const double arg = r.d * D + r.m * M + r.mp * Mp + r.f * F + r.om * Om;
        d_psi += (r.psi_a + r.psi_at * T) * std::sin(arg);
        d_eps += (r.eps_b + r.eps_bt * T) * std::cos(arg);
    }
    return {d_psi * 1e-4, d_eps * 1e-4};
}

}  // namespace spp::model::iau1980
