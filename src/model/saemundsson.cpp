#include "spp/model/saemundsson.hpp"

#include <algorithm>
#include <cmath>

#include "spp/core/constants.hpp"

namespace spp::model::saemundsson {

double apparent_altitude_deg(double h_geo, double pressure_mbar, double temp_c) {
    // Allowed slightly below the horizon for continuity at sunset (~-0.83 deg);
    // a hard cutoff at -0.57 deg would cause a ~0.8 deg step at sunrise.
    if (h_geo < -2.0) return h_geo;

    const double kelvin = 273.0 + temp_c;
    const double scaling = (pressure_mbar / 1010.0) * (283.0 / kelvin);

    // Solve h_app = h_geo + R(h_app), with R = 1.02' / tan(h + 10.3 / (h + 5.11)).
    double h_app = h_geo;
    for (int i = 0; i < 3; ++i) {
        const double h_used = std::max(h_app, -1.0);

        const double div_inner = h_used + 5.11;
        if (div_inner == 0.0) break;

        const double tan_arg_deg = h_used + 10.3 / div_inner;
        const double tan_arg_rad = tan_arg_deg * deg2rad;
        if (std::abs(std::cos(tan_arg_rad)) < 1e-9) break;

        const double r0_arcmin = 1.02 / std::tan(tan_arg_rad);
        const double r_arcmin = r0_arcmin * scaling;
        const double r_deg = r_arcmin / 60.0;

        const double h_new = h_geo + r_deg;
        if (std::abs(h_new - h_app) < 1e-6) {
            h_app = h_new;
            break;
        }
        h_app = h_new;
    }
    return h_app;
}

}  // namespace spp::model::saemundsson
