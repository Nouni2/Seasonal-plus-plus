#include "spp/sun.hpp"

#include "earth_state.hpp"
#include "spp/core/angles.hpp"
#include "spp/sky.hpp"

namespace spp::sun {

Ecliptic geometric(Instant<TT> t) {
    // The geocentric Sun is the heliocentric Earth (VSOP87D) reflected through
    // the origin: λ = L + 180°, β = -B, same distance.
    const auto e = detail::earth_state(t);
    return {norm360(e.l_rad * rad2deg + 180.0), -(e.b_rad * rad2deg), e.r_au};
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

Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm) {
    const auto tt = t.tt();
    const auto n = earth::nutation(tt);
    return sky::observe(apparent(tt, n), obs, earth::lst_deg(t.ut1(), n, obs.lon_deg), atm);
}

}  // namespace spp::sun
