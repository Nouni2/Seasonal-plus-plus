// pipeline_check: fidelity of the whole pipeline against an extended-precision
// evaluation of the same models.
//
// The reference re-implements every step in x87 long double (64-bit mantissa)
// from the models' own coefficients, without calling the library's code: time
// scales, VSOP87D (complete series), IAU 1980 nutation, Laskar obliquity, the
// FK5 correction, aberration, the ecliptic-to-equatorial transform, GMST/GAST,
// WGS84 parallax, diurnal aberration, horizontal coordinates, Saemundsson
// refraction, and the events (transit, sunrise, sunset, civil twilight) solved
// to 10 ns. Two inputs come from the library: TAI - UTC (a table lookup) and
// ΔT outside the leap-second era (a polynomial whose double rounding is below
// 1e-11 s). It checks that the code computes its models faithfully, including
// the Chebyshev cache and every reformulation; it does not check that the
// models are right (tools/reference_check.py does).
//
// Reports max, RMS and bias over 1000 BC-3000 AD for RA, Dec, altitude,
// azimuth and event instants, prints PASS or FAIL, and exits with 0 or 1.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <print>
#include <random>
#include <span>
#include <string>
#include <vector>

#include "model/iau1980_terms.hpp"
#include "spp/model/espenak_meeus.hpp"
#include "spp/model/vsop87d.hpp"
#include "spp/spp.hpp"
#include "spp/vsop87d_earth_terms.inc"

namespace {

using ld = long double;
namespace vsop = spp::model::vsop87d;

static_assert(std::numeric_limits<ld>::digits > std::numeric_limits<double>::digits,
              "pipeline_check needs an extended-precision long double");

constexpr ld kPi = 3.141592653589793238462643383279502884L;
constexpr ld kDeg = kPi / 180.0L;
constexpr ld kArcsec = 3600.0L;

ld norm360(ld x) {
    const ld r = fmodl(x, 360.0L);
    return r < 0 ? r + 360.0L : r;
}

ld wrap180(ld x) { return norm360(x + 180.0L) - 180.0L; }

// ---------------------------------------------------------------------------
// Time: an instant as days since J2000.0 on each scale, in long double.

struct Scales {
    ld utc, ut1, tt;   // days since JD 2451545.0
};

Scales scales(std::int64_t utc_day, ld utc_frac, double dut1_s) {
    const ld utc = ld(utc_day - spp::jd_j2000_day) + utc_frac;
    // TAI - UTC is constant over the UTC day; look it up for the day of this instant.
    const ld whole = floorl(utc_frac);
    const double frac_in_day = double(utc_frac - whole);
    const auto tai = spp::tai_minus_utc({utc_day + std::int64_t(whole), frac_in_day});
    ld tt_minus_utc;
    if (tai) {
        tt_minus_utc = 32.184L + ld(*tai);
    } else {
        const double year = 2000.0 + double(utc / 365.25L);
        tt_minus_utc = ld(dut1_s) + ld(spp::model::espenak_meeus::delta_t_seconds(year));
    }
    return {utc, utc + ld(dut1_s) / 86400.0L, utc + tt_minus_utc / 86400.0L};
}

// ---------------------------------------------------------------------------
// Earth state and the apparent Sun.

ld series(const std::span<const vsop::Term> (&s)[6], ld t) {
    ld value = 0.0L, t_power = 1.0L;
    for (const auto& terms : s) {
        ld inner = 0.0L;
        for (const auto& term : terms) inner += ld(term.a) * cosl(ld(term.b) + ld(term.c) * t);
        value += inner * t_power;
        t_power *= t;
    }
    return value;
}

struct Nutation {
    ld dpsi_deg, deps_deg;
};

Nutation nutation(ld T) {
    const ld T2 = T * T, T3 = T2 * T;
    const ld D = norm360(ld(297.85036) + ld(445267.111480) * T - ld(0.0019142) * T2 + T3 / ld(189474.0)) * kDeg;
    const ld M = norm360(ld(357.52772) + ld(35999.050340) * T - ld(0.0001603) * T2 - T3 / ld(300000.0)) * kDeg;
    const ld Mp = norm360(ld(134.96298) + ld(477198.867398) * T + ld(0.0086972) * T2 + T3 / ld(56250.0)) * kDeg;
    const ld F = norm360(ld(93.27191) + ld(483202.017538) * T - ld(0.0036825) * T2 + T3 / ld(327270.0)) * kDeg;
    const ld Om = norm360(ld(125.04452) - ld(1934.136261) * T + ld(0.0020708) * T2 + T3 / ld(450000.0)) * kDeg;
    ld dpsi = 0.0L, deps = 0.0L;
    for (const auto& r : spp::model::iau1980::kTerms) {
        const ld arg = r.d * D + r.m * M + r.mp * Mp + r.f * F + r.om * Om;
        dpsi += (ld(r.psi_a) + ld(r.psi_at) * T) * sinl(arg);
        deps += (ld(r.eps_b) + ld(r.eps_bt) * T) * cosl(arg);
    }
    return {dpsi * 1e-4L / kArcsec, deps * 1e-4L / kArcsec};
}

ld mean_obliquity_deg(ld T) {
    const ld U = T / 100.0L;
    const ld c[] = {ld(84381.448), ld(-4680.93), ld(-1.55), ld(1999.25), ld(-51.38), ld(-249.67),
                    ld(-39.05),    ld(7.12),     ld(27.87), ld(5.79),    ld(2.45)};
    ld v = 0.0L;
    for (int i = 10; i >= 0; --i) v = v * U + c[i];
    return v / kArcsec;
}

struct Sun {
    ld ra_deg, dec_deg, r_au;
    ld dpsi_deg, eps_deg;
};

Sun apparent_sun(ld tt_days) {
    const ld T = tt_days / 36525.0L;
    const ld t = T / 10.0L;
    const ld L = series(vsop::data::L, t), B = series(vsop::data::B, t), R = series(vsop::data::R, t);
    const Nutation n = nutation(T);
    const ld eps = mean_obliquity_deg(T) + n.deps_deg;

    // Geocentric = heliocentric Earth reflected; VSOP87 dynamical frame -> FK5.
    const ld lon0 = norm360(L / kDeg + 180.0L);
    const ld lat0 = -B / kDeg;
    const ld lp = (lon0 - ld(1.397) * T - ld(0.00031) * T * T) * kDeg;
    const ld lon = lon0 - ld(0.09033) / kArcsec;
    const ld lat = lat0 + ld(0.03916) / kArcsec * (cosl(lp) - sinl(lp));

    // Nutation in longitude and annual aberration.
    const ld lambda = (lon + n.dpsi_deg - ld(spp::sun_aberration_arcsec) / R / kArcsec) * kDeg;
    const ld beta = lat * kDeg, e = eps * kDeg;
    const ld dec = asinl(sinl(beta) * cosl(e) + cosl(beta) * sinl(e) * sinl(lambda));
    ld ra = atan2l(sinl(lambda) * cosl(e) - tanl(beta) * sinl(e), cosl(lambda));
    if (ra < 0) ra += 2.0L * kPi;
    return {ra / kDeg, dec / kDeg, R, n.dpsi_deg, eps};
}

// Local apparent sidereal time, degrees.
ld lst_deg(ld ut1_days, const Sun& s, double lon_deg) {
    const ld T = ut1_days / 36525.0L;
    const ld gmst = ld(280.46061837) + 360.0L * ut1_days + ld(0.98564736629) * ut1_days + ld(0.000387933) * T * T
                    + ld(-1.0 / 38710000.0) * T * T * T;
    return norm360(gmst + ld(lon_deg) + s.dpsi_deg * cosl(s.eps_deg * kDeg));
}

// ---------------------------------------------------------------------------
// The observer.

struct Topo {
    ld alt_geom_deg, az_deg, alt_refr_deg, upper_limb_deg;
};

ld refraction(ld h, ld pressure_mbar) {
    const ld scaling = (pressure_mbar / 1010.0L) * (283.0L / (273.0L + 10.0L));
    auto r = [](ld x) { return ld(1.02) / tanl((x + ld(10.3) / (x + ld(5.11))) * kDeg) + ld(0.0019279); };
    if (h >= -2.0L) return h + r(h) * scaling / 60.0L;
    if (h <= -4.0L) return h;
    return h + r(-2.0L) * ((h + 4.0L) / 2.0L) * scaling / 60.0L;
}

Topo observe(const Sun& s, ld lst, const spp::Observer& obs) {
    const ld phi = ld(obs.lat_deg) * kDeg;
    const ld f = ld(spp::wgs84_f), a = ld(spp::wgs84_a_m);
    const ld u = atanl((1.0L - f) * tanl(phi));
    const ld rc = cosl(u) + ld(obs.elev_m) / a * cosl(phi);
    const ld rs = (1.0L - f) * sinl(u) + ld(obs.elev_m) / a * sinl(phi);

    // Parallax.
    const ld sin_pi = sinl(ld(spp::sun_parallax_arcsec) / s.r_au / kArcsec * kDeg);
    const ld alpha = s.ra_deg * kDeg, delta = s.dec_deg * kDeg;
    const ld H = lst * kDeg - alpha;
    const ld den = cosl(delta) - rc * sin_pi * cosl(H);
    const ld d_alpha = atan2l(-rc * sin_pi * sinl(H), den);
    const ld delta1 = atan2l((sinl(delta) - rs * sin_pi) * cosl(d_alpha), den);
    const ld alpha1 = alpha + d_alpha;

    // Diurnal aberration.
    const ld k = ld(spp::earth_rotation_rad_s) * a / ld(spp::speed_of_light_m_s) * rc;
    const ld H1 = lst * kDeg - alpha1;
    const ld alpha2 = alpha1 + k * cosl(H1) / cosl(delta1);
    const ld delta2 = delta1 + k * sinl(H1) * sinl(delta1);

    // Horizontal.
    const ld H2 = wrap180((lst * kDeg - alpha2) / kDeg) * kDeg;
    const ld alt = asinl(sinl(phi) * sinl(delta2) + cosl(phi) * cosl(delta2) * cosl(H2));
    ld az = atan2l(-cosl(delta2) * sinl(H2), sinl(delta2) * cosl(phi) - cosl(delta2) * cosl(H2) * sinl(phi));
    if (az < 0) az += 2.0L * kPi;

    const ld h = std::clamp(ld(obs.elev_m), -1000.0L, 11000.0L);
    const ld pressure = 1010.0L * powl(1.0L - ld(2.25577e-5) * h, ld(5.25588));
    const ld alt_deg = alt / kDeg;
    const ld sd = ld(spp::sun_radius_arcsec) / s.r_au / kArcsec;
    return {alt_deg, az / kDeg, refraction(alt_deg, pressure), refraction(alt_deg + sd, pressure)};
}

struct State {
    Sun sun;
    ld lst;
    Topo topo;
};

State at(std::int64_t utc_day, ld utc_frac, double dut1_s, const spp::Observer& obs) {
    const Scales sc = scales(utc_day, utc_frac, dut1_s);
    const Sun s = apparent_sun(sc.tt);
    const ld lst = lst_deg(sc.ut1, s, obs.lon_deg);
    return {s, lst, observe(s, lst, obs)};
}

// Root of g near x0 (fraction of the UTC day), by the secant method to 1e-13 day.
template <class G>
std::optional<ld> solve(G g, ld x0) {
    ld a = x0, b = x0 + 1e-5L;
    ld fa = g(a), fb = g(b);
    for (int i = 0; i < 60; ++i) {
        if (fb == fa) break;
        const ld c = b - fb * (b - a) / (fb - fa);
        a = b;
        fa = fb;
        b = c;
        fb = g(b);
        if (fabsl(b - a) < 1e-13L) return b;
        if (fabsl(b - x0) > 0.1L) return std::nullopt;
    }
    return fabsl(b - a) < 1e-11L ? std::optional<ld>(b) : std::nullopt;
}

// ---------------------------------------------------------------------------
// Statistics.

struct Stat {
    std::string name;
    std::string unit;
    double limit;
    double max = 0.0, sum2 = 0.0, sum = 0.0;
    long n = 0;

    void add(double e) {
        max = std::max(max, std::fabs(e));
        sum2 += e * e;
        sum += e;
        ++n;
    }
    bool print() const {
        const double rms = n ? std::sqrt(sum2 / double(n)) : 0.0;
        const double bias = n ? sum / double(n) : 0.0;
        const bool ok = max <= limit;
        std::println("{:<40}{:>8}{:>12.3e}{:>12.3e}{:>12.3e}  {:<7}{}", name, n, max, rms, bias, unit,
                     ok ? "" : "FAIL");
        return ok;
    }
};

double seconds(const spp::Time& t, ld ref_frac) {
    return double((ld(t.utc().frac) - ref_frac) * 86400.0L);
}

}  // namespace

int main() {
    std::mt19937_64 rng(1001);
    std::uniform_real_distribution<double> jd_dist(1355818.5, 2816787.5);   // 1000 BC .. 3000 AD
    std::uniform_real_distribution<double> lat_dist(-89.0, 89.0), lon_dist(-180.0, 180.0), elev_dist(0.0, 3000.0);
    std::uniform_real_distribution<double> dut1_dist(-0.9, 0.9);

    // Limits: about twice the rounding level measured for Spp 1.0.1 (wiki: Accuracy, Numerical precision),
    // so that any change that loses precision fails.
    Stat ra{"right ascension x cos(dec)", "arcsec", 2e-5}, dec{"declination", "arcsec", 2e-5};
    Stat alt{"altitude, refracted", "arcsec", 2e-5}, alt_g{"altitude, geometric", "arcsec", 2e-5};
    Stat az{"azimuth x cos(alt)", "arcsec", 2e-5};
    Stat transit{"transit", "s", 1.5e-6};
    Stat rise{"sunrise, |lat| <= 60", "s", 1.5e-6}, set{"sunset, |lat| <= 60", "s", 1.5e-6};
    Stat rise_hi{"sunrise, |lat| > 60", "s", 1e-5}, set_hi{"sunset, |lat| > 60", "s", 1e-5};
    Stat dusk{"civil dusk, |lat| <= 55", "s", 2e-6};
    long unresolved = 0;

    // Positions at random instants and observers.
    for (int i = 0; i < 20000; ++i) {
        const double jd = jd_dist(rng);
        const auto day = std::int64_t(std::floor(jd));
        const double frac = jd - double(day);
        const spp::Observer obs{lat_dist(rng), lon_dist(rng), elev_dist(rng)};
        const double dut1 = dut1_dist(rng);
        const auto t = spp::Time::from_jd(day, frac, dut1);
        const auto eq = spp::sun::apparent(t.tt());
        const auto hor = spp::sun::horizontal(t, obs);
        const auto hor_g = spp::sun::horizontal(t, obs, {.refraction = false});
        const State ref = at(t.utc().day, ld(t.utc().frac), dut1, obs);
        ra.add(double(wrap180(ld(eq.ra_deg) - ref.sun.ra_deg) * cosl(ref.sun.dec_deg * kDeg) * kArcsec));
        dec.add(double((ld(eq.dec_deg) - ref.sun.dec_deg) * kArcsec));
        alt.add(double((ld(hor.alt_deg) - ref.topo.alt_refr_deg) * kArcsec));
        alt_g.add(double((ld(hor_g.alt_deg) - ref.topo.alt_geom_deg) * kArcsec));
        az.add(double(wrap180(ld(hor.az_deg) - ref.topo.az_deg) * cosl(ref.topo.alt_geom_deg * kDeg) * kArcsec));
    }

    // Events on random dates and observers.
    std::uniform_int_distribution<int> year_dist(-1000, 2999), month_dist(1, 12), mday_dist(1, 28);
    for (int i = 0; i < 4000; ++i) {
        const spp::CivilDate date{year_dist(rng), month_dist(rng), mday_dist(rng)};
        if (date.year == 1582 && date.month == 10) continue;
        const spp::Observer obs{lat_dist(rng), lon_dist(rng), elev_dist(rng)};
        const double dut1 = dut1_dist(rng);
        const auto d = spp::events::day(date, obs, {}, dut1);

        const auto tr = d.transit.utc();
        const auto hour_angle = [&](ld x) {
            const State s = at(tr.day, x, dut1, obs);
            return wrap180(s.lst - s.sun.ra_deg);
        };
        if (const auto x = solve(hour_angle, ld(tr.frac))) transit.add(seconds(d.transit, *x));
        else ++unresolved;

        const bool high = std::fabs(obs.lat_deg) > 60.0;
        for (const auto& [ev, stat] : {std::pair{d.rise, high ? &rise_hi : &rise}, std::pair{d.set, high ? &set_hi : &set}}) {
            if (!ev) continue;
            const auto u = ev->utc();
            const auto limb = [&](ld x) { return at(u.day, x, dut1, obs).topo.upper_limb_deg; };
            if (const auto x = solve(limb, ld(u.frac))) stat->add(seconds(*ev, *x));
            else ++unresolved;
        }

        if (std::fabs(obs.lat_deg) <= 55.0) {
            const auto c = spp::events::crossing(d.transit, obs, {}, spp::events::civil_twilight_deg,
                                                 spp::events::Direction::After);
            if (c) {
                const auto u = c->utc();
                const auto civil = [&](ld x) { return at(u.day, x, dut1, obs).topo.alt_refr_deg + 6.0L; };
                if (const auto x = solve(civil, ld(u.frac))) dusk.add(seconds(*c, *x));
                else ++unresolved;
            }
        }
    }

    std::println("Spp minus an extended-precision evaluation of the same models, 1000 BC .. 3000 AD");
    std::println("{:<40}{:>8}{:>12}{:>12}{:>12}", "", "n", "max", "RMS", "bias");
    bool pass = true;
    for (const Stat* s : {&ra, &dec, &alt_g, &alt, &az, &transit, &rise, &set, &rise_hi, &set_hi, &dusk})
        pass = s->print() && pass;
    std::println("reference roots not resolved: {}", unresolved);
    std::println("{}", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
