#include "spp/sun.hpp"

#include <cmath>

#include "build_mode.hpp"
#include "earth_state.hpp"
#include "spp/core/angles.hpp"
#include "spp/sky.hpp"

namespace spp::sun {

Ecliptic geometric(Instant<TT> t) {
    // 1. Earth heliocentric (VSOP87D), degrees with longitude in [0, 360).
    const auto e = detail::earth_state(t);

    if constexpr (!detail::kExact) {
        // The geocentric Sun is the heliocentric Earth reflected through the
        // origin: λ = L + 180°, β = -B, same distance. (Going through
        // rectangular coordinates and back, as below, only adds rounding.)
        return {norm360(e.l_rad * rad2deg + 180.0), -(e.b_rad * rad2deg), e.r_au};
    }

    const double L_E_rad = norm360(e.l_rad * rad2deg) * deg2rad;
    const double B_E_rad = (e.b_rad * rad2deg) * deg2rad;
    const double R_E_au = e.r_au;

    // 2-3. Heliocentric Earth vector, inverted to the geocentric Sun.
    const double X_S = -(R_E_au * std::cos(B_E_rad) * std::cos(L_E_rad));
    const double Y_S = -(R_E_au * std::cos(B_E_rad) * std::sin(L_E_rad));
    const double Z_S = -(R_E_au * std::sin(B_E_rad));

    const double R_S_au = std::sqrt(X_S * X_S + Y_S * Y_S + Z_S * Z_S);

    double lambda_rad = std::atan2(Y_S, X_S);
    const double beta_rad = std::atan2(Z_S, std::hypot(X_S, Y_S));
    if (lambda_rad < 0.0) lambda_rad += two_pi;

    return {lambda_rad * rad2deg, beta_rad * rad2deg, R_S_au};
}

Ecliptic apparent_ecliptic(Instant<TT> t, const earth::Nutation& n) {
    const Ecliptic geo = geometric(t);

    // Nutation in longitude, then annual aberration -20.4898" / R (Meeus 25.10).
    // Aberration in latitude is < 0.01" for the Sun and is neglected.
    const double aberration_deg = (sun_aberration_arcsec / geo.dist_au) / 3600.0;
    const double lambda_deg = norm360(geo.lon_deg + n.dpsi_deg - aberration_deg);
    return {lambda_deg, geo.lat_deg, geo.dist_au};
}

Equatorial apparent(Instant<TT> t, const earth::Nutation& n) {
    return sky::to_equatorial(apparent_ecliptic(t, n), n.true_obliquity_deg);
}

Equatorial apparent(Instant<TT> t) { return apparent(t, earth::nutation(t)); }

}  // namespace spp::sun
