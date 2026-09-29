#include "spp/earth.hpp"

#include <cmath>

#include "earth_state.hpp"
#include "spp/core/angles.hpp"
#include "spp/model/laskar.hpp"
#include "spp/model/meeus.hpp"

namespace spp::earth {

Nutation nutation(Instant<TT> t) {
    const auto s = detail::earth_state(t);

    const double dpsi_deg = s.dpsi_arcsec / 3600.0;
    const double deps_deg = s.deps_arcsec / 3600.0;
    const double eps0_deg = model::laskar::mean_obliquity_deg(t.centuries_since_j2000());
    return {dpsi_deg, deps_deg, eps0_deg, eps0_deg + deps_deg};
}

double mean_obliquity_deg(Instant<TT> t) {
    return model::laskar::mean_obliquity_deg(t.centuries_since_j2000());
}

double delta_t(const Time& t) { return t.delta_t(); }

double gmst_deg(Instant<UT1> t) { return model::meeus::gmst_deg(t.day - jd_j2000_day, t.frac); }

double equation_of_equinoxes_deg(const Nutation& n) {
    return n.dpsi_deg * std::cos(n.true_obliquity_deg * deg2rad);
}

double gast_deg(Instant<UT1> t, const Nutation& n) {
    return norm360(gmst_deg(t) + equation_of_equinoxes_deg(n));
}

double gast_deg(const Time& t) { return gast_deg(t.ut1(), nutation(t.tt())); }

double lst_deg(Instant<UT1> t, const Nutation& n, double lon_east_deg) {
    return norm360(gmst_deg(t) + lon_east_deg + equation_of_equinoxes_deg(n));
}

double lst_deg(const Time& t, double lon_east_deg) {
    return lst_deg(t.ut1(), nutation(t.tt()), lon_east_deg);
}

}  // namespace spp::earth
