#pragma once

// ΔT = TT - UT1 in seconds: observed monthly values from IERS EOP 20 C04,
// 1972-01 to 2026-09 (linearly interpolated), a bridge to 2050, and Espenak &
// Meeus (NASA/GSFC) polynomials elsewhere.
namespace spp::model::espenak_meeus {

double delta_t_seconds(double decimal_year);

}  // namespace spp::model::espenak_meeus
