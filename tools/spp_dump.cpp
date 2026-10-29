// spp_dump: prints Spp results for requests read from standard input, one
// output line per request line, every number with 17 significant digits. It is
// the library side of tools/reference_check.py, which compares the results
// with independent references.
//
// Requests and output columns ("nan" stands for an unset value):
//
//   pos <jd_day> <utc_frac> <dut1_s|nan> <lat_deg> <lon_deg> <elev_m>
//       tai_utc tt_utc dut1 delta_t tt_day tt_frac ut1_day ut1_frac gmst gast
//       dpsi deps eps0 eps geo_lon geo_lat geo_r app_lon ra dec alt_geom az
//       alt_refr semidiameter                                   (angles in degrees)
//   tt  <jd_day> <tt_frac>              the instant is passed as given, unnormalized
//       ra dec geo_lon geo_r app_lon app_lat
//   day <year> <month> <day> <lat_deg> <lon_deg> <elev_m> <dut1_s|nan>
//       type transit_day transit_frac rise_day rise_frac set_day set_frac duration_h
//   alt <jd_day> <utc_frac> <dut1_s|nan> <lat_deg> <lon_deg> <elev_m>
//       alt_geom semidiameter alt_refr upper_limb_alt
//   cross <year> <month> <day> <lat_deg> <lon_deg> <elev_m> <dut1_s|nan> <target_deg> <-1|1>
//       crossing_day crossing_frac               (center at target, before/after transit)
//   cal <year> <month> <day> <hour> <minute> <second>
//       jd_day jd_frac year month day hour minute second iso_utc
//   dt  <decimal_year>
//       delta_t_s                                     (the ΔT model, spp::model::espenak_meeus)
//
// Instants in the output are UTC split Julian Dates unless named otherwise.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>

#include "spp/model/espenak_meeus.hpp"
#include "spp/spp.hpp"

namespace {

std::optional<double> opt(double v) { return std::isnan(v) ? std::nullopt : std::optional<double>(v); }

// Reads one number; accepts "nan", which operator>> rejects.
double num(std::istringstream& in) {
    std::string s;
    in >> s;
    return std::stod(s);
}

long long whole(std::istringstream& in) { return static_cast<long long>(num(in)); }

void put(double v) { std::printf(" %.17g", v); }

void put_time(const std::optional<spp::Time>& t) {
    if (!t) {
        std::printf(" nan nan");
        return;
    }
    const auto u = t->utc();
    std::printf(" %lld %.17g", static_cast<long long>(u.day), u.frac);
}

void pos(std::istringstream& in) {
    const long long day = whole(in);
    const double frac = num(in), dut1 = num(in), lat = num(in), lon = num(in), elev = num(in);
    const auto t = spp::Time::from_jd(day, frac, opt(dut1));
    const auto tt = t.tt();
    const auto ut1 = t.ut1();
    const auto n = spp::earth::nutation(tt);
    const auto geo = spp::sun::geometric(tt);
    const auto app = spp::sun::apparent_ecliptic(tt, n);
    const auto eq = spp::sun::apparent(tt, n);
    const spp::Observer obs{lat, lon, elev};
    const auto hor_geom = spp::sun::horizontal(t, obs, spp::Atmosphere{.refraction = false});
    const auto hor_refr = spp::sun::horizontal(t, obs);
    const auto tai = spp::tai_minus_utc(t.utc());

    std::printf("pos");
    put(tai.value_or(std::numeric_limits<double>::quiet_NaN()));
    put(t.tt_minus_utc());
    put(t.dut1());
    put(t.delta_t());
    std::printf(" %lld", static_cast<long long>(tt.day));
    put(tt.frac);
    std::printf(" %lld", static_cast<long long>(ut1.day));
    put(ut1.frac);
    put(spp::earth::gmst_deg(ut1));
    put(spp::earth::gast_deg(ut1, n));
    put(n.dpsi_deg);
    put(n.deps_deg);
    put(n.mean_obliquity_deg);
    put(n.true_obliquity_deg);
    put(geo.lon_deg);
    put(geo.lat_deg);
    put(geo.dist_au);
    put(app.lon_deg);
    put(eq.ra_deg);
    put(eq.dec_deg);
    put(hor_geom.alt_deg);
    put(hor_geom.az_deg);
    put(hor_refr.alt_deg);
    put(spp::sun::semidiameter_deg(tt));
    std::printf("\n");
}

void tt_state(std::istringstream& in) {
    const long long day = whole(in);
    const double frac = num(in);
    const spp::Instant<spp::TT> tt{day, frac};
    const auto n = spp::earth::nutation(tt);
    const auto eq = spp::sun::apparent(tt, n);
    const auto geo = spp::sun::geometric(tt);
    const auto ecl = spp::sun::apparent_ecliptic(tt, n);
    std::printf("tt");
    put(eq.ra_deg);
    put(eq.dec_deg);
    put(geo.lon_deg);
    put(geo.dist_au);
    put(ecl.lon_deg);
    put(ecl.lat_deg);
    std::printf("\n");
}

void day(std::istringstream& in) {
    const int y = static_cast<int>(whole(in)), m = static_cast<int>(whole(in)), d = static_cast<int>(whole(in));
    const double lat = num(in), lon = num(in), elev = num(in), dut1 = num(in);
    const spp::Observer obs{lat, lon, elev};
    const auto r = spp::events::day({y, m, d}, obs, {}, opt(dut1));
    std::printf("day %d", static_cast<int>(r.type));
    put_time(r.transit);
    put_time(r.rise);
    put_time(r.set);
    put(r.duration_hours);
    std::printf("\n");
}

void alt(std::istringstream& in) {
    const long long day = whole(in);
    const double frac = num(in), dut1 = num(in), lat = num(in), lon = num(in), elev = num(in);
    const auto t = spp::Time::from_jd(day, frac, opt(dut1));
    const spp::Observer obs{lat, lon, elev};
    std::printf("alt");
    put(spp::sun::horizontal(t, obs, spp::Atmosphere{.refraction = false}).alt_deg);
    put(spp::sun::semidiameter_deg(t.tt()));
    put(spp::sun::horizontal(t, obs).alt_deg);
    put(spp::events::upper_limb_altitude_deg(t, obs, {}));
    std::printf("\n");
}

void cross(std::istringstream& in) {
    const int y = static_cast<int>(whole(in)), m = static_cast<int>(whole(in)), d = static_cast<int>(whole(in));
    const double lat = num(in), lon = num(in), elev = num(in), dut1 = num(in), target = num(in);
    const auto dir = num(in) < 0.0 ? spp::events::Direction::Before : spp::events::Direction::After;
    const spp::Observer obs{lat, lon, elev};
    const auto r = spp::events::day({y, m, d}, obs, {}, opt(dut1));
    std::printf("cross");
    put_time(spp::events::crossing(r.transit, obs, {}, target, dir));
    std::printf("\n");
}

void cal(std::istringstream& in) {
    const int y = static_cast<int>(whole(in)), mo = static_cast<int>(whole(in)), d = static_cast<int>(whole(in));
    const int h = static_cast<int>(whole(in)), mi = static_cast<int>(whole(in));
    const double s = num(in);
    const auto t = spp::Time::from_utc(y, mo, d, h, mi, s, 0.0);
    const auto u = t.utc();
    const auto c = t.to_civil();
    std::printf("cal %lld %.17g %d %d %d %d %d %.17g %s\n", static_cast<long long>(u.day), u.frac, c.year, c.month,
                c.day, c.hour, c.minute, c.second, t.iso_utc().c_str());
}

void dt(std::istringstream& in) {
    std::printf("dt");
    put(spp::model::espenak_meeus::delta_t_seconds(num(in)));
    std::printf("\n");
}

}  // namespace

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string kind;
        in >> kind;
        if (kind == "pos") pos(in);
        else if (kind == "tt") tt_state(in);
        else if (kind == "day") day(in);
        else if (kind == "alt") alt(in);
        else if (kind == "cross") cross(in);
        else if (kind == "cal") cal(in);
        else if (kind == "dt") dt(in);
        else std::printf("?\n");
    }
}
