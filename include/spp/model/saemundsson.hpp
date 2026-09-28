#pragma once

// Saemundsson (1986) refraction (Meeus 16.4), scaled for pressure and
// temperature and solved iteratively for the apparent altitude.
namespace spp::model::saemundsson {

double apparent_altitude_deg(double geometric_alt_deg, double pressure_mbar, double temp_c);

}  // namespace spp::model::saemundsson
