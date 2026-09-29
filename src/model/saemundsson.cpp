#include "spp/model/saemundsson.hpp"

#include <cmath>

#include "spp/core/constants.hpp"

namespace spp::model::saemundsson {

namespace {

// R = 1.02' / tan(h + 10.3° / (h + 5.11°)) + 0.0019279', h the true (geometric)
// altitude in degrees. The constant makes R = 0 at the zenith (Meeus 16.4).
double refraction_arcmin(double h) { return 1.02 / std::tan((h + 10.3 / (h + 5.11)) * deg2rad) + 0.0019279; }

// The formula applies down to kFormulaFloor. Below it the refraction fades
// linearly to 0 at kFadeEnd, so the apparent altitude stays continuous and
// increasing, and twilight altitudes (-6° and lower) remain geometric.
constexpr double kFormulaFloor = -2.0;
constexpr double kFadeEnd = -4.0;

}  // namespace

double apparent_altitude_deg(double h_geo, double pressure_mbar, double temp_c) {
    const double scaling = (pressure_mbar / 1010.0) * (283.0 / (273.0 + temp_c));
    if (h_geo >= kFormulaFloor) return h_geo + refraction_arcmin(h_geo) * scaling / 60.0;
    if (h_geo <= kFadeEnd) return h_geo;
    const double weight = (h_geo - kFadeEnd) / (kFormulaFloor - kFadeEnd);
    return h_geo + refraction_arcmin(kFormulaFloor) * weight * scaling / 60.0;
}

}  // namespace spp::model::saemundsson
