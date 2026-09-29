#include "spp/sky.hpp"

#include <algorithm>
#include <cmath>

#include "spp/core/angles.hpp"
#include "spp/core/constants.hpp"
#include "spp/model/saemundsson.hpp"

namespace spp::sky {

Equatorial to_equatorial(const Ecliptic& ecl, double true_obliquity_deg) {
    const double eps_rad = true_obliquity_deg * deg2rad;
    const double lon_rad = ecl.lon_deg * deg2rad;
    const double lat_rad = ecl.lat_deg * deg2rad;

    const double sin_lon = std::sin(lon_rad);
    const double cos_lon = std::cos(lon_rad);
    const double sin_lat = std::sin(lat_rad);
    const double cos_lat = std::cos(lat_rad);
    const double sin_eps = std::sin(eps_rad);
    const double cos_eps = std::cos(eps_rad);

    // sin δ = sin β cos ε + cos β sin ε sin λ
    const double sin_delta = std::clamp((sin_lat * cos_eps) + (cos_lat * sin_eps * sin_lon), -1.0, 1.0);
    const double delta_rad = std::asin(sin_delta);

    // tan α = (sin λ cos ε - tan β sin ε) / cos λ
    const double y = (sin_lon * cos_eps) - (std::tan(lat_rad) * sin_eps);
    double alpha_rad = std::atan2(y, cos_lon);
    if (alpha_rad < 0.0) alpha_rad += two_pi;

    return {alpha_rad * rad2deg, delta_rad * rad2deg, ecl.dist_au};
}

Horizontal to_horizontal(const Equatorial& eq, double lat_deg, double lst_deg) {
    const double lat_rad = lat_deg * deg2rad;
    const double dec_rad = eq.dec_deg * deg2rad;
    const double ha_rad = wrap180(lst_deg - eq.ra_deg) * deg2rad;   // H = LST - α

    const double sin_lat = std::sin(lat_rad);
    const double cos_lat = std::cos(lat_rad);
    const double sin_dec = std::sin(dec_rad);
    const double cos_dec = std::cos(dec_rad);
    const double sin_h = std::sin(ha_rad);
    const double cos_h = std::cos(ha_rad);

    // sin h = sin φ sin δ + cos φ cos δ cos H
    const double sin_alt = std::clamp((sin_lat * sin_dec) + (cos_lat * cos_dec * cos_h), -1.0, 1.0);
    const double alt_rad = std::asin(sin_alt);

    // Azimuth from North, increasing East.
    const double y = -cos_dec * sin_h;
    const double x = (sin_dec * cos_lat) - (cos_dec * cos_h * sin_lat);
    double az_rad = std::atan2(y, x);
    if (az_rad < 0.0) az_rad += two_pi;

    return {alt_rad * rad2deg, az_rad * rad2deg};
}

namespace {

// Observer's geocentric position in units of the equatorial radius: distance
// from the axis (ρ cos φ') and from the equatorial plane (ρ sin φ').
struct GeocentricObserver {
    double rho_cos_phi;
    double rho_sin_phi;
};

GeocentricObserver geocentric_observer(const Observer& obs) {
    const double lat_rad = obs.lat_deg * deg2rad;
    const double u = std::atan((1.0 - wgs84_f) * std::tan(lat_rad));   // reduced latitude
    const double h_ratio = obs.elev_m / wgs84_a_m;
    return {std::cos(u) + h_ratio * std::cos(lat_rad), (1.0 - wgs84_f) * std::sin(u) + h_ratio * std::sin(lat_rad)};
}

}  // namespace

Equatorial parallax(const Equatorial& geocentric, const Observer& obs, double lst_deg) {
    const auto [rho_cos_phi, rho_sin_phi] = geocentric_observer(obs);

    // Horizontal parallax at the Sun's distance.
    const double pi_rad = ((sun_parallax_arcsec / geocentric.dist_au) / 3600.0) * deg2rad;
    const double sin_pi = std::sin(pi_rad);

    const double ra_rad = geocentric.ra_deg * deg2rad;
    const double h_rad = lst_deg * deg2rad - ra_rad;
    const double sin_h = std::sin(h_rad);
    const double cos_h = std::cos(h_rad);

    const double dec_rad = geocentric.dec_deg * deg2rad;
    const double cos_dec = std::cos(dec_rad);
    const double sin_dec = std::sin(dec_rad);

    // Meeus 40.2 / 40.3, with atan2 for both.
    const double numerator = -rho_cos_phi * sin_pi * sin_h;
    const double denominator = cos_dec - (rho_cos_phi * sin_pi * cos_h);
    const double delta_alpha_rad = std::atan2(numerator, denominator);

    const double num_dec = (sin_dec - rho_sin_phi * sin_pi) * std::cos(delta_alpha_rad);
    const double dec_topo_rad = std::atan2(num_dec, denominator);

    double ra_topo_rad = ra_rad + delta_alpha_rad;
    if (ra_topo_rad < 0.0) ra_topo_rad += two_pi;
    else if (ra_topo_rad >= two_pi) ra_topo_rad -= two_pi;

    return {ra_topo_rad * rad2deg, dec_topo_rad * rad2deg, geocentric.dist_au};
}

Equatorial diurnal_aberration(const Equatorial& topocentric, const Observer& obs, double lst_deg) {
    // The observer moves eastward at ω a ρ cos φ'. To first order in v/c the
    // apparent direction shifts toward the east point:
    //   Δα = κ ρ cos φ' cos H / cos δ,  Δδ = κ ρ cos φ' sin H sin δ,
    // with κ = ω a / c ≈ 0.320".
    constexpr double kappa = earth_rotation_rad_s * wgs84_a_m / speed_of_light_m_s;
    const double k = kappa * geocentric_observer(obs).rho_cos_phi;
    const double h_rad = (lst_deg - topocentric.ra_deg) * deg2rad;
    const double dec_rad = topocentric.dec_deg * deg2rad;
    const double d_alpha = k * std::cos(h_rad) / std::cos(dec_rad);
    const double d_delta = k * std::sin(h_rad) * std::sin(dec_rad);
    return {norm360(topocentric.ra_deg + d_alpha * rad2deg), topocentric.dec_deg + d_delta * rad2deg,
            topocentric.dist_au};
}

Atmosphere at_elevation(const Atmosphere& atm, double elev_m) {
    if (atm.pressure_mbar) return atm;
    // International Standard Atmosphere pressure ratio (troposphere), applied
    // to the 1010 mbar standard that the refraction model assumes.
    const double h = std::clamp(elev_m, -1000.0, 11000.0);
    Atmosphere scaled = atm;
    scaled.pressure_mbar = 1010.0 * std::pow(1.0 - 2.25577e-5 * h, 5.25588);
    return scaled;
}

double refraction(double geometric_alt_deg, const Atmosphere& atm) {
    if (!atm.refraction) return geometric_alt_deg;
    return model::saemundsson::apparent_altitude_deg(geometric_alt_deg, atm.pressure_mbar.value_or(1010.0),
                                                     atm.temp_c);
}

Horizontal observe(const Equatorial& geocentric, const Observer& obs, double lst_deg, const Atmosphere& atm) {
    const Equatorial topo = diurnal_aberration(parallax(geocentric, obs, lst_deg), obs, lst_deg);
    const Horizontal hor = to_horizontal(topo, obs.lat_deg, lst_deg);
    return {refraction(hor.alt_deg, at_elevation(atm, obs.elev_m)), hor.az_deg};
}

}  // namespace spp::sky
