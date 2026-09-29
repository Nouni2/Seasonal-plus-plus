// Throughput of the solar pipeline: single instants and a date x latitude grid.
// Every measurement runs once untimed first, so all of them are timed with the
// ephemeris cache already filled for the instants they use.

#include <chrono>
#include <print>

#include "spp/spp.hpp"

namespace {

using Clock = std::chrono::steady_clock;

volatile double sink = 0.0;

// Runs `body` once to warm up, then again timed; returns seconds.
template <class F>
double timed(F body) {
    body();
    const auto t0 = Clock::now();
    body();
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

}  // namespace

int main() {
    constexpr int n = 20000;
    const spp::Observer paris{48.8566, 2.3522, 35.0};
    const auto instant = [](int i) { return spp::Time::from_jd(2461000, i * 0.0137); };

    // 1. Apparent geocentric Sun (VSOP87D + nutation + aberration).
    const double s1 = timed([&] {
        for (int i = 0; i < n; ++i) sink = sink + spp::sun::apparent(instant(i).tt()).ra_deg;
    });
    std::println("sun::apparent          {:>10.2f} us/call  ({} calls)", s1 / n * 1e6, n);

    // 2. Full topocentric refracted altitude (includes the apparent Sun).
    const double s2 = timed([&] {
        for (int i = 0; i < n; ++i) sink = sink + spp::events::altitude_deg(instant(i), paris, {});
    });
    std::println("events::altitude_deg   {:>10.2f} us/call  ({} calls)", s2 / n * 1e6, n);

    // 3. Sunrise/sunset for every day of a year at 19 latitudes.
    int days = 0;
    const double s3 = timed([&] {
        days = 0;
        for (int lat = -90; lat <= 90; lat += 10) {
            auto t = spp::Time::from_utc(2026, 1, 1);
            for (int d = 0; d < 365; ++d, ++days) {
                const auto c = t.to_civil();
                sink = sink + spp::events::day({c.year, c.month, c.day}, {double(lat), 0.0}, {}).duration_hours;
                t = t.shifted(1.0);
            }
        }
    });
    std::println("events::day            {:>10.2f} us/day   ({} days, {:.3f} s total)", s3 / days * 1e6, days, s3);
    return 0;
}
