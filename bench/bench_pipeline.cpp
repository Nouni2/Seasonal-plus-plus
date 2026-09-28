// Throughput of the solar pipeline: single instants and a date x latitude grid.

#include <chrono>
#include <print>

#include "spp/spp.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

}  // namespace

int main() {
    volatile double sink = 0.0;

    // 1. Apparent geocentric Sun (VSOP87D + nutation + aberration).
    {
        constexpr int n = 20000;
        const auto t0 = Clock::now();
        for (int i = 0; i < n; ++i) {
            const auto t = spp::Time::from_jd(2461000, i * 0.0137);
            sink = sink + spp::sun::apparent(t.tt()).ra_deg;
        }
        const double s = seconds_since(t0);
        std::println("sun::apparent          {:>10.2f} us/call  ({} calls)", s / n * 1e6, n);
    }

    // 2. Full topocentric refracted altitude.
    {
        constexpr int n = 20000;
        const spp::Observer paris{48.8566, 2.3522, 35.0};
        const auto t0 = Clock::now();
        for (int i = 0; i < n; ++i) {
            const auto t = spp::Time::from_jd(2461000, i * 0.0137);
            sink = sink + spp::events::altitude_deg(t, paris, {});
        }
        const double s = seconds_since(t0);
        std::println("events::altitude_deg   {:>10.2f} us/call  ({} calls)", s / n * 1e6, n);
    }

    // 3. Sunrise/sunset for every day of a year at 19 latitudes.
    {
        int days = 0;
        const auto t0 = Clock::now();
        for (int lat = -90; lat <= 90; lat += 10) {
            auto t = spp::Time::from_utc(2026, 1, 1);
            for (int d = 0; d < 365; ++d, ++days) {
                const auto c = t.to_civil();
                const auto r = spp::events::day({c.year, c.month, c.day}, {double(lat), 0.0}, {});
                sink = sink + r.duration_hours;
                t = t.shifted(1.0);
            }
        }
        const double s = seconds_since(t0);
        std::println("events::day            {:>10.2f} us/day   ({} days, {:.3f} s total)", s / days * 1e6,
                     days, s);
    }
    return 0;
}
