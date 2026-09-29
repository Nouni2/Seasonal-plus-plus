#pragma once

#include "spp/core/types.hpp"

// The observer's sky: frame transforms and observer-specific corrections.
namespace spp::sky {

Equatorial to_equatorial(const Ecliptic& ecl, double true_obliquity_deg);

// Geometric transform only: no parallax, no refraction.
Horizontal to_horizontal(const Equatorial& eq, double lat_deg, double lst_deg);

// Geocentric -> topocentric equatorial (WGS84 observer position).
Equatorial parallax(const Equatorial& geocentric, const Observer& obs, double lst_deg);

// Diurnal aberration: the shift from the observer's rotational velocity
// (up to 0.32"), applied to a topocentric position.
Equatorial diurnal_aberration(const Equatorial& topocentric, const Observer& obs, double lst_deg);

// `atm` with its pressure filled in: when unset, the standard 1010 mbar scaled
// to elevation `elev_m` with the International Standard Atmosphere.
Atmosphere at_elevation(const Atmosphere& atm, double elev_m);

// Geometric altitude -> apparent (refracted) altitude. An unset pressure is
// taken as 1010 mbar (sea level).
double refraction(double geometric_alt_deg, const Atmosphere& atm);

// What the observer sees: parallax, diurnal aberration, horizontal
// coordinates, then refraction of the altitude with the pressure at the
// observer's elevation. `geocentric` is an apparent geocentric position.
Horizontal observe(const Equatorial& geocentric, const Observer& obs, double lst_deg, const Atmosphere& atm);

}  // namespace spp::sky
