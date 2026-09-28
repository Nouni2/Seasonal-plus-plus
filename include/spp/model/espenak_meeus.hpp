#pragma once

// ΔT = TT - UT1 in seconds: IERS/USNO observed values 2005-2026 (linearly
// interpolated), a bridge to 2050, and Espenak & Meeus (NASA/GSFC)
// polynomials elsewhere.
namespace spp::model::espenak_meeus {

double delta_t_seconds(double decimal_year);

}  // namespace spp::model::espenak_meeus
