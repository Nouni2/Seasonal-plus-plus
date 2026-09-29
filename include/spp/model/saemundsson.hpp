#pragma once

// Saemundsson (1986) refraction (Meeus 16.4), a function of the true
// (geometric) altitude, scaled for pressure and temperature. Applied down to
// -2°, then faded linearly to 0 at -4° so the apparent altitude is continuous.
namespace spp::model::saemundsson {

double apparent_altitude_deg(double geometric_alt_deg, double pressure_mbar, double temp_c);

}  // namespace spp::model::saemundsson
