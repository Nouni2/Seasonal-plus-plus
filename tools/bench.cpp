// Throughput of the solar pipeline.
//
// Warm path: consecutive instants and a date x latitude grid. Each of these
// measurements runs once untimed first, so it is timed with the ephemeris cache
// already filled for the instants it uses.
// Cold path: random instants and dates over 1000 BC .. 3000 AD, timed on first
// use, so almost every call builds a new 16-day cache segment.

#include <chrono>
#include <cmath>
#include <print>
#include <random>
#include <vector>

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

// Times `body` on its first run; returns seconds.
template <class F>
double timed_cold(F body) {
    const auto t0 = Clock::now();
    body();
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

#ifdef NDEBUG
constexpr const char* kBuild = "release (NDEBUG)";
#else
constexpr const char* kBuild = "debug";
#endif

}  // namespace

int main() {
    std::println("compiler {}, {} build; one thread", __VERSION__, kBuild);
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

    // 4. Cold: the full topocentric altitude at random instants, 1000 BC .. 3000 AD.
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> jd(1355818.5, 2816787.5);
    constexpr int n_cold = 2000;
    std::vector<spp::Time> cold;
    for (int i = 0; i < n_cold; ++i) {
        const double x = jd(rng);
        cold.push_back(spp::Time::from_jd(static_cast<std::int64_t>(std::floor(x)), x - std::floor(x)));
    }
    const double s4 = timed_cold([&] {
        for (const auto& t : cold) sink = sink + spp::events::altitude_deg(t, paris, {});
    });
    std::println("events::altitude_deg   {:>10.2f} us/call  cold, random instants ({} calls)", s4 / n_cold * 1e6,
                 n_cold);

    // 5. Cold: events::day on random dates, 1000 BC .. 3000 AD.
    std::uniform_int_distribution<int> year(-1000, 2999), month(1, 12), mday(1, 28);
    std::vector<spp::CivilDate> dates;
    while (static_cast<int>(dates.size()) < n_cold) {
        const spp::CivilDate d{year(rng), month(rng), mday(rng)};
        if (spp::is_valid(d)) dates.push_back(d);
    }
    const double s5 = timed_cold([&] {
        for (const auto& d : dates) sink = sink + spp::events::day(d, paris).duration_hours;
    });
    std::println("events::day            {:>10.2f} us/day   cold, random dates ({} days)", s5 / n_cold * 1e6, n_cold);
    return 0;
}
