#pragma once

#include <numbers>

namespace spp {

inline constexpr double pi = std::numbers::pi;
inline constexpr double two_pi = 2.0 * std::numbers::pi;

inline constexpr double deg2rad = pi / 180.0;
inline constexpr double rad2deg = 180.0 / pi;
inline constexpr double hours2deg = 15.0;

// Epochs and time units
inline constexpr long long jd_j2000_day = 2451545;   // JD 2451545.0 = 2000-01-01 12:00 TT
inline constexpr double days_per_century = 36525.0;
inline constexpr double seconds_per_day = 86400.0;
inline constexpr long long jd_gregorian_start = 2299161;   // 1582-10-15

// WGS84 ellipsoid and Earth rotation
inline constexpr double wgs84_a_m = 6378137.0;
inline constexpr double wgs84_f = 1.0 / 298.257223563;
inline constexpr double earth_rotation_rad_s = 7.292115e-5;   // sidereal rate
inline constexpr double speed_of_light_m_s = 299792458.0;

// Solar constants
inline constexpr double sun_aberration_arcsec = 20.4898;       // Meeus 25.10, divided by R
inline constexpr double sun_parallax_arcsec = 8.794;           // at 1 AU
inline constexpr double sun_radius_arcsec = 959.63;            // semidiameter at 1 AU
inline constexpr double sun_semidiameter_deg = 16.0 / 60.0;    // mean, almanac convention

}  // namespace spp
