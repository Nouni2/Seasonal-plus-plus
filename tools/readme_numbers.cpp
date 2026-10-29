// readme_numbers: computes, with the library, the figures the README quotes
// about events, sidereal time and time scales, so each can be reproduced.
// Figures measured against external references are produced by
// tools/reference_check.py; accuracy and speed figures by spp_cache_check,
// spp_pipeline_check and spp_bench.

#include <cmath>
#include <print>

#include "spp/model/espenak_meeus.hpp"
#include "spp/spp.hpp"

namespace {

using spp::events::Direction;

double seconds_between(const spp::Time& a, const spp::Time& b) {
    const auto x = a.utc();
    const auto y = b.utc();
    return (static_cast<double>(x.day - y.day) + (x.frac - y.frac)) * 86400.0;
}

std::string hhmm(const spp::Time& t) {
    const auto c = t.to_civil();
    return std::format("{:02d}:{:02d}:{:02.0f}", c.hour, c.minute, std::floor(c.second));
}

}  // namespace

int main() {
    const spp::Observer paris{48.8566, 2.3522, 35.0};

    // Gallery: day length on 20 March 2026 (figure observers: sea level, longitude 0 for the equator).
    std::println("[gallery] day length 2026-03-20: equator {:.2f} h, Paris {:.2f} h",
                 spp::events::day({2026, 3, 20}, {0.0, 0.0, 0.0}).duration_hours,
                 spp::events::day({2026, 3, 20}, paris).duration_hours);

    // Gallery: solar noon in Paris at the June solstice, March equinox and December solstice.
    std::println("[gallery] Paris solar noon (UTC): {}, {}, {}",
                 hhmm(spp::events::day({2026, 6, 21}, paris).transit),
                 hhmm(spp::events::day({2026, 3, 20}, paris).transit),
                 hhmm(spp::events::day({2026, 12, 21}, paris).transit));

    // Gallery: azimuth range of sunrise and sunset in Paris over 2026.
    double rise_lo = 400, rise_hi = -1, set_lo = 400, set_hi = -1;
    for (auto t = spp::Time::from_utc(2026, 1, 1, 12); t.to_civil().year == 2026; t = t.shifted(1.0)) {
        const auto c = t.to_civil();
        const auto d = spp::events::day({c.year, c.month, c.day}, paris);
        if (d.rise) {
            const double az = spp::sun::horizontal(*d.rise, paris).az_deg;
            rise_lo = std::fmin(rise_lo, az);
            rise_hi = std::fmax(rise_hi, az);
        }
        if (d.set) {
            const double az = spp::sun::horizontal(*d.set, paris).az_deg;
            set_lo = std::fmin(set_lo, az);
            set_hi = std::fmax(set_hi, az);
        }
    }
    std::println("[gallery] Paris 2026 azimuth: sunrise {:.1f}..{:.1f} deg, sunset {:.1f}..{:.1f} deg", rise_lo,
                 rise_hi, set_lo, set_hi);

    // §13.1: the upper limb condition in geometric terms, default atmosphere, 16' semidiameter.
    double lo = -2.0, hi = 0.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        (spp::sky::refraction(mid, {}) < 0.0 ? lo : hi) = mid;
    }
    const double limb = 0.5 * (lo + hi), center = limb - 16.0 / 60.0;
    std::println("[13.1] upper limb at geometric {:.4f} deg (refraction {:.2f}'), center {:.4f} deg "
                 "(refraction there {:.2f}')",
                 limb, -limb * 60.0, center, (spp::sky::refraction(center, {}) - center) * 60.0);

    // §13.1: offsets from the almanac convention (geometric center at -0.8333 deg), sea level, 2026.
    const struct {
        const char* name;
        spp::CivilDate date;
        spp::Observer obs;
    } rows[] = {{"equator, March equinox", {2026, 3, 20}, {0.0, 0.0, 0.0}},
                {"48.86 N (Paris), June solstice", {2026, 6, 21}, {48.8566, 2.3522, 0.0}},
                {"48.86 N (Paris), December solstice", {2026, 12, 21}, {48.8566, 2.3522, 0.0}},
                {"60 N, June solstice", {2026, 6, 21}, {60.0, 0.0, 0.0}}};
    for (const auto& r : rows) {
        const auto d = spp::events::day(r.date, r.obs);
        const spp::Atmosphere none{.refraction = false};
        const auto a = spp::events::crossing(d.transit, r.obs, none, -50.0 / 60.0, Direction::Before);
        const auto b = spp::events::crossing(d.transit, r.obs, none, -50.0 / 60.0, Direction::After);
        std::println("[13.1] almanac offset, {}: sunrise {:+.1f} s, sunset {:+.1f} s", r.name,
                     seconds_between(*d.rise, *a), seconds_between(*d.set, *b));
    }

    // §13.2: a day next to polar day.
    {
        const auto d = spp::events::day({2026, 7, 27}, {70.0, 0.0, 0.0});
        std::println("[13.2] 70 N 2026-07-27: sunrise {}, sunset {}, duration {:.2f} h",
                     d.rise ? d.rise->iso_utc() : "-", d.set ? d.set->iso_utc() : "-", d.duration_hours);
    }

    // §13.3: rate of the Sun's geocentric hour angle over 2026, degrees per day.
    double rate_lo = 1e9, rate_hi = -1e9;
    for (auto t = spp::Time::from_utc(2026, 1, 1); t.to_civil().year == 2026; t = t.shifted(0.25)) {
        const auto h = [](const spp::Time& u) {
            return spp::earth::gast_deg(u) - spp::sun::apparent(u.tt()).ra_deg;
        };
        const double dt = 1.0 / 24.0;
        const double rate = std::remainder(h(t.shifted(dt)) - h(t), 360.0) / dt;
        rate_lo = std::fmin(rate_lo, rate);
        rate_hi = std::fmax(rate_hi, rate);
    }
    std::println("[13.3] hour angle rate over 2026: {:.3f} .. {:.3f} deg/day", rate_lo, rate_hi);

    // §2.5: TT - UTC across 1972-01-01 0h UTC (model before, leap seconds after).
    const auto before = spp::Time::from_jd(2441317, 0.5 - 1e-9);
    const auto after = spp::Time::from_jd(2441317, 0.5);
    std::println("[2.5] TT - UTC at 1972-01-01: {:.4f} s -> {:.4f} s, step {:.4f} s", before.tt_minus_utc(),
                 after.tt_minus_utc(), before.tt_minus_utc() - after.tt_minus_utc());

    // §2.5: the DUT1 held after the leap-second table.
    std::println("[2.5] DUT1 held after 2027-06-28: {:.4f} s", spp::Time::from_utc(2028, 1, 1).dut1());

    // §2.6: the ΔT bridge.
    std::println("[2.6] ΔT bridge: {:.2f} s in 2030, {:.2f} s in 2040",
                 spp::model::espenak_meeus::delta_t_seconds(2030.0), spp::model::espenak_meeus::delta_t_seconds(2040.0));
}
