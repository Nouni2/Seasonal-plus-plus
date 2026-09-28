#pragma once

#include "spp/core/types.hpp"

// The observer's sky: frame transforms and observer-specific corrections.
namespace spp::sky {

Equatorial to_equatorial(const Ecliptic& ecl, double true_obliquity_deg);

// Geometric transform only: no parallax, no refraction.
Horizontal to_horizontal(const Equatorial& eq, double lat_deg, double lst_deg);

// Geocentric -> topocentric equatorial (WGS84 observer position).
Equatorial parallax(const Equatorial& geocentric, const Observer& obs, double lst_deg);

// Geometric altitude -> apparent (refracted) altitude.
double refraction(double geometric_alt_deg, const Atmosphere& atm);

}  // namespace spp::sky
