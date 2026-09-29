// Solar events: transit by Newton steps on the hour angle; crossings from an
// analytic first guess, a tight bracket and Illinois refinement to
// SolverConfig::time_tolerance_s (1 µs).

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
    return sky::observe(g.sun, obs, earth::lst_deg(g.ut1, g.nutation, obs.lon_deg), atm).alt_deg;
}

double hour_angle_deg(const Geocentric& g, double lon_east_deg) {
    return wrap180(earth::lst_deg(g.ut1, g.nutation, lon_east_deg) - g.sun.ra_deg);
}

Time approx_noon(const CivilDate& date, double lon_east_deg, std::optional<double> dut1_s) {
    // 12:00 - lon/15 in UTC, split into hours, minutes and seconds.
    const double approx_noon_utc_hour = 12.0 - (lon_east_deg / 15.0);
    const int h = static_cast<int>(approx_noon_utc_hour);
    const double rem_minutes = (approx_noon_utc_hour - h) * 60.0;
    const int m = static_cast<int>(rem_minutes);
    const double s = (rem_minutes - m) * 60.0;
    return Time::from_utc(date.year, date.month, date.day, h, m, s, dut1_s);
}

// Time above the horizon within the 24 h centered on transit. A missing rise
// (or set) means the Sun is already up at the start (or still up at the end)
// of that window, which happens next to polar day.
double duration_hours(const Time& transit, const std::optional<Time>& rise, const std::optional<Time>& set) {
    const auto start = (rise ? *rise : transit.shifted(-0.5)).utc();
    const auto end = (set ? *set : transit.shifted(0.5)).utc();
    return (static_cast<double>(end.day - start.day) + (end.frac - start.frac)) * 24.0;
}

struct Transit {
    Time time;
    Geocentric state;   // at `time`
};

Transit transit_impl(const Time& approx, double lon_east_deg, const SolverConfig& cfg) {
    // Newton on the hour angle: it advances ~360° per day, so each step
    // shrinks the error ~3000x. 1e-9° of hour angle is ~2e-7 s.
    Time current = approx;
    Geocentric g = geocentric(current);
    double H = hour_angle_deg(g, lon_east_deg);
    for (int i = 0; i < cfg.max_iter && std::abs(H) >= 1e-9; ++i) {
        const Time next = current.shifted(-H / 360.0);
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
    // Instants on the transit day, as offsets from transit (keeps DUT1 handling).
    auto at = [&](double frac) { return tr.time.shifted(frac - base.frac); };
    auto f = [&](double frac) { return altitude_deg(at(frac), obs, atm) - target_alt_deg; };

    const double sgn = static_cast<int>(dir);
    double lo = base.frac, f_lo = alt_transit - target_alt_deg;   // f >= 0 side
    if (f_lo < 0.0) return std::nullopt;                          // transit below target
    double hi = 0.0, f_hi = 0.0;                                  // f < 0 side
    bool bracketed = false;

    // 1. First guess from the hour angle where the geometric altitude reaches
    //    the target, with refraction at the target folded in:
    //    cos H0 = (sin h0 - sin φ sin δ) / (cos φ cos δ)
    const double refr = sky::refraction(target_alt_deg, sky::at_elevation(atm, obs.elev_m)) - target_alt_deg;
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

}  // namespace

double altitude_deg(const Time& t, const Observer& obs, const Atmosphere& atm) {
    return altitude_from(geocentric(t), obs, atm);
}

Time transit(const Time& approx, double lon_east_deg, const SolverConfig& cfg) {
    return transit_impl(approx, lon_east_deg, cfg).time;
}

std::optional<Time> crossing(const Time& transit_time, const Observer& obs, const Atmosphere& atm,
                             double target_alt_deg, Direction dir, const SolverConfig& cfg) {
    const Geocentric g = geocentric(transit_time);
    return crossing_impl({transit_time, g}, altitude_from(g, obs, atm), obs, atm, target_alt_deg, dir, cfg);
}

Day day(const CivilDate& date, const Observer& obs, const Atmosphere& atm, std::optional<double> dut1_s,
        const SolverConfig& cfg) {
    const Transit tr = transit_impl(approx_noon(date, obs.lon_deg, dut1_s), obs.lon_deg, cfg);

    // The upper limb touches the horizon when the center is one semidiameter
    // (959.63" / R, taken at transit) below it.
    const double limb = -sun::semidiameter_deg(tr.time.tt());

    // Classify from the altitude at transit and at both lower culminations
    // around it (approximated as 12 h before and after).
    const double alt_noon = altitude_from(tr.state, obs, atm);
    if (alt_noon < limb) return {date, DayType::PolarNight, tr.time, {}, {}, 0.0};
    const double alt_before = altitude_deg(tr.time.shifted(-0.5), obs, atm);
    const double alt_after = altitude_deg(tr.time.shifted(0.5), obs, atm);
    if (alt_before > limb && alt_after > limb) return {date, DayType::PolarDay, tr.time, {}, {}, 24.0};

    auto rise = crossing_impl(tr, alt_noon, obs, atm, limb, Direction::Before, cfg);
    auto set = crossing_impl(tr, alt_noon, obs, atm, limb, Direction::After, cfg);
    return {date, DayType::Normal, tr.time, rise, set, duration_hours(tr.time, rise, set)};
}

}  // namespace spp::events
