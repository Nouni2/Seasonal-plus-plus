#include "observe.hpp"

#include "spp/sky.hpp"
#include "spp/sun.hpp"

namespace spp::detail {

Geocentric geocentric(const Time& t) {
    const auto tt = t.tt();
    const auto n = earth::nutation(tt);
    return {n, sun::apparent(tt, n), t.ut1()};
}

double lst_deg(const Geocentric& g, double lon_east_deg) { return earth::lst_deg(g.ut1, g.nutation, lon_east_deg); }

Horizontal observed(const Geocentric& g, const Observer& obs, const Atmosphere& atm) {
    return sky::observe(g.sun, obs, lst_deg(g, obs.lon_deg), atm);
}

}  // namespace spp::detail
