// Two solvers share the altitude pipeline:
//   - default: analytic first guess, tight bracket, Illinois refinement to
//     SolverConfig::time_tolerance_s (1 µs), transit refined to ~1e-7 s;
//   - SPP_EXACT (reference mode): hourly scan from transit, then
//     secant/bisection to 0.001 deg of altitude.

#include "spp/events.hpp"

#include <cmath>

#include "spp/core/angles.hpp"
#include "spp/earth.hpp"
#include "spp/sky.hpp"
#include "spp/sun.hpp"

namespace spp::events {

namespace {

// Everything about one instant that does not depend on the observer.
// Nutation is computed once and shared by the apparent place and the
// equation of the equinoxes, so both always use the same Δψ.
struct Geocentric {
    earth::Nutation nutation;
    Equatorial sun;
    Instant<UT1> ut1;
};

Geocentric geocentric(const Time& t) {
    const auto tt = t.tt();
    const auto n = earth::nutation(tt);
    return {n, sun::apparent(tt, n), t.ut1()};
}

double altitude_from(const Geocentric& g, const Observer& obs, const Atmosphere& atm) {
    const double lst = earth::lst_deg(g.ut1, g.nutation, obs.lon_deg);
    const Equatorial topo = sky::parallax(g.sun, obs, lst);
    const Horizontal hor = sky::to_horizontal(topo, obs.lat_deg, lst);
    return sky::refraction(hor.alt_deg, atm);
}

double hour_angle_deg(const Geocentric& g, double lon_east_deg) {
    return wrap180(earth::lst_deg(g.ut1, g.nutation, lon_east_deg) - g.sun.ra_deg);
}

Time approx_noon(const CivilDate& date, double lon_east_deg, double dut1_s) {
    // 12:00 - lon/15 in UTC, split into hours, minutes and seconds.
    const double approx_noon_utc_hour = 12.0 - (lon_east_deg / 15.0);
    const int h = static_cast<int>(approx_noon_utc_hour);
    const double rem_minutes = (approx_noon_utc_hour - h) * 60.0;
    const int m = static_cast<int>(rem_minutes);
    const double s = (rem_minutes - m) * 60.0;
    return Time::from_utc(date.year, date.month, date.day, h, m, s, dut1_s);
}

double duration_hours(const std::optional<Time>& rise, const std::optional<Time>& set) {
    if (!rise || !set) return 0.0;
    const auto r = rise->utc();
    const auto st = set->utc();
    return (static_cast<double>(st.day - r.day) + (st.frac - r.frac)) * 24.0;
}

#ifdef SPP_EXACT

// ---------------------------------------------------------------------------
// Reference-mode solver
// ---------------------------------------------------------------------------

constexpr double kEpsilonDeg = 0.001;
constexpr int kMaxIter = 30;

Time transit_impl(const Time& approx, double lon_east_deg) {
    Time current = approx;
    for (int i = 0; i < 3; ++i) {
        const double H = hour_angle_deg(geocentric(current), lon_east_deg);
        if (std::abs(H) < 0.001) break;
        current = Time::from_jd(current.utc().day, current.utc().frac - H / 360.0, current.dut1());
    }
    return current;
}

std::optional<Time> crossing_impl(const Time& transit_time, const Observer& obs, const Atmosphere& atm,
                                  double target_alt_deg, Direction dir) {
    const auto base = transit_time.utc();
    const double dut1 = transit_time.dut1();
    auto at = [&](double frac) { return Time::from_jd(base.day, frac, dut1); };
    auto f = [&](double frac) { return altitude_deg(at(frac), obs, atm) - target_alt_deg; };

    // 1. Bracket by scanning outward from transit in 1 h steps (max 12 h).
    double frac_above = base.frac;
    double val_above = f(frac_above);
    if (val_above < 0.0) return std::nullopt;

    std::optional<double> frac_below;
    double val_below = 0.0;
    const double step = static_cast<int>(dir) * (1.0 / 24.0);
    for (int i = 1; i <= 12; ++i) {
        const double frac_i = base.frac + i * step;
        const double val_i = f(frac_i);
        if (val_i < 0.0) {
            frac_below = frac_i;
            val_below = val_i;
            break;
        }
        frac_above = frac_i;
        val_above = val_i;
    }
    if (!frac_below) return std::nullopt;

    // 2. Secant/bisection hybrid; the bracket is kept at every step.
    double lo = frac_above, f_lo = val_above;
    double hi = *frac_below, f_hi = val_below;
    for (int iter = 0; iter < kMaxIter; ++iter) {
        const double denom = f_hi - f_lo;
        double cand = std::abs(denom) > 1e-12 ? hi - f_hi * (hi - lo) / denom : 0.5 * (lo + hi);

        const double lo_b = lo < hi ? lo : hi;
        const double hi_b = lo < hi ? hi : lo;
        if (!(lo_b < cand && cand < hi_b)) cand = 0.5 * (lo + hi);

        const double f_cand = f(cand);
        if (std::abs(f_cand) < kEpsilonDeg) return at(cand);

        if (f_cand > 0.0) {
            lo = cand;
            f_lo = f_cand;
        } else {
            hi = cand;
            f_hi = f_cand;
        }
        if (std::abs(hi - lo) < 1.0 / 86400.0 / 10.0) return at(0.5 * (lo + hi));
    }
    return std::nullopt;
}

#else

// ---------------------------------------------------------------------------
// Default solver
// ---------------------------------------------------------------------------

struct Transit {
    Time time;
    Geocentric state;   // at `time`
};

Transit transit_impl(const Time& approx, double lon_east_deg, const SolverConfig& cfg) {
    // Newton on the hour angle: it advances ~360° per day, so each step
    // shrinks the error ~1000x. 1e-9° of hour angle is ~2e-7 s.
    Time current = approx;
    Geocentric g = geocentric(current);
    double H = hour_angle_deg(g, lon_east_deg);
    for (int i = 0; i < cfg.max_iter && std::abs(H) >= 1e-9; ++i) {
        const Time next = Time::from_jd(current.utc().day, current.utc().frac - H / 360.0, current.dut1());
        const Geocentric g_next = geocentric(next);
        const double H_next = hour_angle_deg(g_next, lon_east_deg);
        if (std::abs(H_next) >= std::abs(H)) break;   // at the noise floor: keep the best
        current = next;
        g = g_next;
        H = H_next;
    }
    return {current, g};
}

std::optional<Time> crossing_impl(const Transit& tr, double alt_transit, const Observer& obs,
                                  const Atmosphere& atm, double target_alt_deg, Direction dir,
                                  const SolverConfig& cfg) {
    const auto base = tr.time.utc();
    const double dut1 = tr.time.dut1();
    auto at = [&](double frac) { return Time::from_jd(base.day, frac, dut1); };
    auto f = [&](double frac) { return altitude_deg(at(frac), obs, atm) - target_alt_deg; };

    const double sgn = static_cast<int>(dir);
    double lo = base.frac, f_lo = alt_transit - target_alt_deg;   // f >= 0 side
    if (f_lo < 0.0) return std::nullopt;                          // transit below target
    double hi = 0.0, f_hi = 0.0;                                  // f < 0 side
    bool bracketed = false;

    // 1. First guess from the hour angle where the geometric altitude reaches
    //    the target, with refraction at the target folded in:
    //    cos H0 = (sin h0 - sin φ sin δ) / (cos φ cos δ)
    const double refr = sky::refraction(target_alt_deg, atm) - target_alt_deg;
    const double h0 = (target_alt_deg - refr) * deg2rad;
    const double phi = obs.lat_deg * deg2rad;
    const double dec = tr.state.sun.dec_deg * deg2rad;
    const double c = (std::sin(h0) - std::sin(phi) * std::sin(dec)) / (std::cos(phi) * std::cos(dec));

    if (std::abs(c) < 1.0) {
        const double guess = base.frac + sgn * (std::acos(c) * rad2deg / 360.0);
        const double fg = f(guess);
        double step = sgn * (2.0 / 1440.0);   // 2 minutes, doubled while walking
        if (fg >= 0.0) {
            lo = guess;
            f_lo = fg;
            while (std::abs(lo + step - base.frac) <= 0.5) {   // walk away from transit
                const double p = lo + step;
                const double fp = f(p);
                if (fp < 0.0) {
                    hi = p;
                    f_hi = fp;
                    bracketed = true;
                    break;
                }
                lo = p;
                f_lo = fp;
                step *= 2.0;
            }
        } else {
            hi = guess;
            f_hi = fg;
            bracketed = true;
            while (true) {   // walk back toward transit, where f >= 0 is known
                const double p = hi - step;
                if (sgn * (p - base.frac) <= 0.0) break;   // reached transit: [transit, hi]
                const double fp = f(p);
                if (fp >= 0.0) {
                    lo = p;
                    f_lo = fp;
                    break;
                }
                hi = p;
                f_hi = fp;
                step *= 2.0;
            }
        }
    }

    // 2. No usable guess (near the polar limits): hourly scan from transit.
    if (!bracketed) {
        lo = base.frac;
        f_lo = alt_transit - target_alt_deg;
        for (int i = 1; i <= 12 && !bracketed; ++i) {
            const double p = base.frac + i * sgn / 24.0;
            const double fp = f(p);
            if (fp < 0.0) {
                hi = p;
                f_hi = fp;
                bracketed = true;
            } else {
                lo = p;
                f_lo = fp;
            }
        }
        if (!bracketed) return std::nullopt;
    }

    // 3. Illinois (modified regula falsi): keeps the bracket and converges
    //    superlinearly. Stops once the bracket is below the time tolerance.
    const double tol_days = cfg.time_tolerance_s / seconds_per_day;
    double x = 0.5 * (lo + hi);
    int side = 0;
    for (int iter = 0; iter < cfg.max_iter && std::abs(hi - lo) > tol_days; ++iter) {
        x = (lo * f_hi - hi * f_lo) / (f_hi - f_lo);
        if (!((lo < x && x < hi) || (hi < x && x < lo))) x = 0.5 * (lo + hi);
        const double fx = f(x);
        if (fx == 0.0) return at(x);
        if (fx > 0.0) {
            lo = x;
            f_lo = fx;
            if (side == 1) f_hi *= 0.5;
            side = 1;
        } else {
            hi = x;
            f_hi = fx;
            if (side == -1) f_lo *= 0.5;
            side = -1;
        }
    }
    return at(0.5 * (lo + hi));
}

#endif

}  // namespace

double altitude_deg(const Time& t, const Observer& obs, const Atmosphere& atm) {
    return altitude_from(geocentric(t), obs, atm);
}

#ifdef SPP_EXACT

Time transit(const Time& approx, double lon_east_deg, const SolverConfig&) {
    return transit_impl(approx, lon_east_deg);
}

std::optional<Time> crossing(const Time& transit_time, const Observer& obs, const Atmosphere& atm,
                             double target_alt_deg, Direction dir, const SolverConfig&) {
    return crossing_impl(transit_time, obs, atm, target_alt_deg, dir);
}

Day day(const CivilDate& date, const Observer& obs, const Atmosphere& atm, double dut1_s,
        const SolverConfig&) {
    const Time t_transit = transit_impl(approx_noon(date, obs.lon_deg, dut1_s), obs.lon_deg);

    const double alt_noon = altitude_deg(t_transit, obs, atm);
    const double alt_midnight = altitude_deg(t_transit.shifted(0.5), obs, atm);
    if (alt_noon < upper_limb_deg) return {date, DayType::PolarNight, t_transit, {}, {}, 0.0};
    if (alt_midnight > upper_limb_deg) return {date, DayType::PolarDay, t_transit, {}, {}, 24.0};

    auto rise = crossing_impl(t_transit, obs, atm, upper_limb_deg, Direction::Before);
    auto set = crossing_impl(t_transit, obs, atm, upper_limb_deg, Direction::After);
    return {date, DayType::Normal, t_transit, rise, set, duration_hours(rise, set)};
}

#else

Time transit(const Time& approx, double lon_east_deg, const SolverConfig& cfg) {
    return transit_impl(approx, lon_east_deg, cfg).time;
}

std::optional<Time> crossing(const Time& transit_time, const Observer& obs, const Atmosphere& atm,
                             double target_alt_deg, Direction dir, const SolverConfig& cfg) {
    const Geocentric g = geocentric(transit_time);
    return crossing_impl({transit_time, g}, altitude_from(g, obs, atm), obs, atm, target_alt_deg, dir, cfg);
}

Day day(const CivilDate& date, const Observer& obs, const Atmosphere& atm, double dut1_s,
        const SolverConfig& cfg) {
    const Transit tr = transit_impl(approx_noon(date, obs.lon_deg, dut1_s), obs.lon_deg, cfg);

    // Classify from the altitude at transit and 12 h later.
    const double alt_noon = altitude_from(tr.state, obs, atm);
    const double alt_midnight = altitude_deg(tr.time.shifted(0.5), obs, atm);
    if (alt_noon < upper_limb_deg) return {date, DayType::PolarNight, tr.time, {}, {}, 0.0};
    if (alt_midnight > upper_limb_deg) return {date, DayType::PolarDay, tr.time, {}, {}, 24.0};

    auto rise = crossing_impl(tr, alt_noon, obs, atm, upper_limb_deg, Direction::Before, cfg);
    auto set = crossing_impl(tr, alt_noon, obs, atm, upper_limb_deg, Direction::After, cfg);
    return {date, DayType::Normal, tr.time, rise, set, duration_hours(rise, set)};
}

#endif

}  // namespace spp::events
