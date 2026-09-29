#include "spp/model/saemundsson.hpp"

#include <cmath>

#include "spp/core/constants.hpp"

namespace spp::model::saemundsson {

double apparent_altitude_deg(double h_geo, double pressure_mbar, double temp_c) {
    // Far below the horizon the formula no longer applies.
    if (h_geo < -2.0) return h_geo;

    // R = 1.02' / tan(h + 10.3° / (h + 5.11°)), h the true (geometric)
    // altitude in degrees, scaled for pressure and temperature.
    const double scaling = (pressure_mbar / 1010.0) * (283.0 / (273.0 + temp_c));
    const double r_arcmin = 1.02 / std::tan((h_geo + 10.3 / (h_geo + 5.11)) * deg2rad);
    return h_geo + r_arcmin * scaling / 60.0;
}

}  // namespace spp::model::saemundsson
