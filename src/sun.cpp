#include "spp/sun.hpp"

#include <cmath>

#include "core/split.hpp"
#include "earth_state.hpp"
#include "observe.hpp"
#include "spp/core/angles.hpp"
#include "spp/sky.hpp"

namespace spp::sun {

Ecliptic geometric(Instant<TT> t) {
    t = detail::normalized(t);
    // The geocentric Sun is the heliocentric Earth (VSOP87D) reflected through
    // the origin: λ = L + 180°, β = -B, same distance.
    const auto e = detail::earth_state(t);
    const double lon = norm360(e.l_rad * rad2deg + 180.0);
    const double lat = -(e.b_rad * rad2deg);

    // From the VSOP87 dynamical ecliptic and equinox to FK5 (Meeus 25.9):
    // λ' = λ - 1.397°T - 0.00031°T², Δλ = -0.09033", Δβ = 0.03916" (cos λ' - sin λ').
    const double T = t.centuries_since_j2000();
    const double lp = (lon - 1.397 * T - 0.00031 * T * T) * deg2rad;
    return {norm360(lon - 0.09033 / 3600.0), lat + 0.03916 / 3600.0 * (std::cos(lp) - std::sin(lp)), e.r_au};
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

double semidiameter_deg(Instant<TT> t) {
    return sun_radius_arcsec / detail::earth_state(detail::normalized(t)).r_au / 3600.0;
}

Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm) {
    return detail::observed(detail::geocentric(t), obs, atm);
}

}  // namespace spp::sun
