#pragma once

// Saemundsson (1986) refraction (Meeus 16.4), a function of the true
// (geometric) altitude, scaled for pressure and temperature.
namespace spp::model::saemundsson {

double apparent_altitude_deg(double geometric_alt_deg, double pressure_mbar, double temp_c);

}  // namespace spp::model::saemundsson
