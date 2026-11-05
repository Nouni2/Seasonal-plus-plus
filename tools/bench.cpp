// Throughput of the solar pipeline and direct/cached lunar calculations.
//
// Solar warm path: consecutive instants and a date x latitude grid. Each of these
// measurements runs once untimed first, so it is timed with the ephemeris cache
// already filled for the instants it uses.
// Solar cold path: random instants and dates over 1000 BC .. 3000 AD, timed on first
// use, so almost every call builds a new 16-day cache segment.

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <format>
#include <print>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "spp/spp.hpp"
#if SPP_BENCH_MOON
#include "spp/model/elp_mpp02.hpp"
#endif

#ifndef SPP_BENCH_CONFIG
#define SPP_BENCH_CONFIG "unspecified"
#endif

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

struct Options {
    bool moon_only = false;
    bool moon_options = false;
    int samples = 128;
    int rounds = 5;
    std::string_view fit = "both";
};

int integer_option(std::string_view text, int low, int high) {
    int value = 0;
    const auto [end, error] = std::from_chars(text.data(),text.data()+text.size(),value);
    if (error != std::errc{} || end != text.data()+text.size() || value < low || value > high)
        throw std::invalid_argument(std::format("benchmark count must be an integer in {}..{}",low,high));
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i=1; i<argc; ++i) {
        const std::string_view arg=argv[i];
        if (arg=="--moon-only") options.moon_only=true;
        else if (arg.starts_with("--samples=")) options.samples=integer_option(arg.substr(10),16,8192);
        else if (arg.starts_with("--rounds=")) options.rounds=integer_option(arg.substr(9),3,31);
        else if (arg.starts_with("--fit=")) {
            options.fit=arg.substr(6);
            if (options.fit!="both" && options.fit!="de405" && options.fit!="llr")
                throw std::invalid_argument("--fit must be de405, llr or both");
        } else throw std::invalid_argument("unknown benchmark option: "+std::string(arg));
        options.moon_options=true;
    }
    return options;
}

const char* compiler_name() {
#if defined(__clang__)
    return "Clang " __clang_version__;
#elif defined(__GNUC__)
    return "GCC " __VERSION__;
#elif defined(_MSC_VER)
    return "MSVC";
#else
    return "unknown compiler";
#endif
}

#if SPP_BENCH_MOON
// Consume every returned numeric field, so LTO cannot turn a full-state or
// full-observation case into a cheaper position/altitude-only calculation.
double checksum(const spp::moon::Position& p) {
    return p.position_km[0]+p.position_km[1]+p.position_km[2];
}
double checksum(const spp::moon::State& s) {
    return s.position_km[0]+s.position_km[1]+s.position_km[2]
        +s.velocity_km_per_day[0]+s.velocity_km_per_day[1]+s.velocity_km_per_day[2];
}
double checksum(const spp::Horizontal& h) { return h.alt_deg+h.az_deg; }
double checksum(const spp::moon::Phase& p) {
    return p.cycle_deg+p.elongation_deg+p.phase_angle_deg+p.illuminated_fraction
        +static_cast<double>(p.waxing);
}
double checksum(const spp::moon::Surface& s) {
    double value=s.subobserver_lon_deg+s.subobserver_lat_deg+s.subsolar_lon_deg+s.subsolar_lat_deg
        +s.pole_position_angle_deg+s.pole_defined;
    for(const auto& row:s.body_to_view)for(double x:row)value+=x;
    return value;
}
double checksum(const spp::moon::Observation& o) {
    return o.equatorial.ra_deg+o.equatorial.dec_deg+o.equatorial.dist_au
        +checksum(o.horizontal)+o.geometric_alt_deg+o.distance_km+o.semidiameter_deg+checksum(o.phase)
        +o.disk.phase_angle_deg+o.disk.illuminated_fraction+o.disk.bright_limb_deg
        +o.disk.parallactic_angle_deg+o.disk.bright_limb_defined+o.disk.zenith_defined+checksum(o.surface);
}

struct Timing {
    double best_us, median_us, worst_us;
};

// A timed sample is an entire batch. File loading, input construction and
// printing stay outside the timer; the checksum is consumed after it stops.
template<class F, class Prepare>
Timing measure(int calls, int rounds, F body, Prepare prepare) {
    prepare();
    sink = sink + body();
    std::vector<double> batches;
    batches.reserve(static_cast<std::size_t>(rounds));
    for (int i=0; i<rounds; ++i) {
        prepare();
        const auto before=Clock::now();
        const double result_checksum=body();
        const double seconds=std::chrono::duration<double>(Clock::now()-before).count();
        sink = sink + result_checksum;
        if (!(seconds>0.0) || !std::isfinite(result_checksum))
            throw std::runtime_error("invalid lunar timing or non-finite result");
        batches.push_back(seconds*1e6/static_cast<double>(calls));
    }
    std::sort(batches.begin(),batches.end());
    const auto middle=batches.size()/2;
    const double median=batches.size()%2 ? batches[middle] : (batches[middle-1]+batches[middle])/2;
    return {batches.front(),median,batches.back()};
}

template<class F>
Timing measure(int calls, int rounds, F body) {
    return measure(calls,rounds,body,[]{});
}

void timing_row(const char* operation, const Timing& t) {
    std::println("  {:<21} {:>11.3f} {:>11.3f} {:>11.3f} {:>12.1f}",
        operation,t.best_us,t.median_us,t.worst_us,1e6/t.median_us);
}

struct LunarInputs {
    const char* name;
    std::vector<spp::Time> utc;
    std::vector<spp::Instant<spp::TDB>> tdb;

    void add(const spp::Time& t) { utc.push_back(t); tdb.push_back(t.tdb()); }
};

void moon_benchmark(const Options& options) {
    LunarInputs nearby{"nearby instants (one-minute steps)",{},{}},
                random{"random dates (2000-2026; seed 7)",{}, {}};
    const auto start=spp::Time::from_utc(2026,11,1,0,0,0.0,0.0);
    const auto search_end=start.shifted(35);
    const auto random_start=spp::Time::from_utc(2000,1,1,0,0,0.0,0.0);
    std::mt19937_64 rng(7);
    for (int i=0; i<options.samples; ++i) {
        nearby.add(start.shifted(static_cast<double>(i)/1440.0));
        const double unit=static_cast<double>(rng()>>11)*0x1p-53;
        random.add(random_start.shifted(unit*9862.0));
    }
    const spp::Observer observer{48.8566,2.3522,35.0};
    std::println("\nMoon: {} queries/batch, {} measured batches, one untimed warmup per case",
                 options.samples,options.rounds);
    std::println("One thread; dates/TDB prepared outside timers; batch means, not latency percentiles.");
    std::println("Direct ELP is compared with twenty-minute Hermite interpolation; both use all coefficients at source evaluations.");
    std::println("Cached-repeat rows warm up the workload. Forced-miss rows clear the lunar cache before EVERY query (clearing included).");
    std::println("Earth/OS caches are not flushed; ELP calls/query exposes capacity misses in repeat workloads.");
    using spp::model::elp_mpp02::Fit;
    for (const Fit fit : {Fit::DE405,Fit::LLR}) {
        const bool de405=fit==Fit::DE405;
        if (options.fit!="both" && options.fit!=(de405?"de405":"llr")) continue;
        const spp::model::elp_mpp02::Ephemeris model(fit);
        spp::model::elp_mpp02::CachedEphemeris cached(fit);
        for (const auto* input : {&nearby,&random}) {
            std::println("\nELP/MPP02 {} | {} | Paris, standard refraction",
                de405?"DE405/DE406":"LLR",input->name);
            std::println("Direct full-series reference:");
            const auto state=measure(options.samples,options.rounds,[&] {
                double sum=0;
                for (const auto& t:input->tdb) sum+=checksum(spp::moon::geometric(t,model));
                return sum;
            });
            const auto position=measure(options.samples,options.rounds,[&] {
                double sum=0;
                for (const auto& t:input->tdb) sum+=checksum(spp::moon::position(t,model));
                return sum;
            });
            const auto horizontal=measure(options.samples,options.rounds,[&] {
                double sum=0;
                for (const auto& t:input->utc) sum+=checksum(spp::moon::horizontal(t,observer,model));
                return sum;
            });
            const auto observe=measure(options.samples,options.rounds,[&] {
                double sum=0;
                for (const auto& t:input->utc) sum+=checksum(spp::moon::observe(t,observer,model));
                return sum;
            });
            const auto phase=measure(options.samples,options.rounds,[&] {
                double sum=0;
                for (const auto& t:input->utc) sum+=checksum(spp::moon::phase(t,model));
                return sum;
            });
            std::println("  {:<21} {:>11} {:>11} {:>11} {:>12}","operation","best us","median us","worst us","queries/s");
            timing_row("geometric state",state);
            timing_row("position only",position);
            timing_row("altitude/azimuth",horizontal);
            timing_row("full observation",observe);
            timing_row("phase/illumination",phase);
            std::println("  Median ratios: state/position {:.3f}x; full observation/coordinates {:.3f}x",
                state.median_us/position.median_us,observe.median_us/horizontal.median_us);

            std::println("Cached provider (same inputs, fit and observation corrections):");
            std::println("  {:<21} {:>11} {:>11} {:>11} {:>12} {:>15}",
                "operation","best us","median us","worst us","queries/s","ELP calls/query");
            auto cached_case=[&](const char* label,const auto& instants,auto operation,bool miss=false) {
                cached.clear();
                double evaluations_per_query=0;
                const auto result=measure(options.samples,options.rounds,[&] {
                    double sum=0;
                    std::uint64_t evaluations=0;
                    for (const auto& t:instants) {
                        if (miss) cached.clear();
                        const auto before=cached.source_evaluations();
                        sum+=checksum(operation(t));
                        evaluations+=cached.source_evaluations()-before;
                    }
                    evaluations_per_query=static_cast<double>(evaluations)/options.samples;
                    return sum;
                });
                std::println("  {:<21} {:>11.3f} {:>11.3f} {:>11.3f} {:>12.1f} {:>15.3f}",label,
                    result.best_us,result.median_us,result.worst_us,1e6/result.median_us,evaluations_per_query);
                return result;
            };
            cached_case("state repeat",input->tdb,[&](auto t){return spp::moon::geometric(t,cached);});
            cached_case("position repeat",input->tdb,[&](auto t){return spp::moon::position(t,cached);});
            const auto cached_horizontal=cached_case("coordinates repeat",input->utc,
                [&](const auto& t){return spp::moon::horizontal(t,observer,cached);});
            cached_case("observation repeat",input->utc,[&](const auto& t){return spp::moon::observe(t,observer,cached);});
            cached_case("phase repeat",input->utc,[&](const auto& t){return spp::moon::phase(t,cached);});
            cached_case("position forced miss",input->tdb,[&](auto t){return spp::moon::position(t,cached);},true);
            const auto missed_horizontal=cached_case("coordinates miss",input->utc,
                [&](const auto& t){return spp::moon::horizontal(t,observer,cached);},true);
            cached_case("observation miss",input->utc,[&](const auto& t){return spp::moon::observe(t,observer,cached);},true);
            std::println("  Coordinates direct/cached: repeat {:.3f}x; forced lunar miss {:.3f}x",
                horizontal.median_us/cached_horizontal.median_us,horizontal.median_us/missed_horizontal.median_us);
        }
        std::size_t events=0;
        const auto search=measure(1,options.rounds,[&] {
            const auto result=spp::moon::phases(start,search_end,model);
            events=result.size();
            double sum=0;
            for (const auto& event:result) sum+=static_cast<double>(event.time.utc().day)
                +event.time.utc().frac+static_cast<int>(event.phase);
            return sum;
        });
        std::println("  Direct 35-day phase search, {} events: best {:.3f}, median {:.3f}, worst {:.3f} ms/search",
            events,search.best_us/1000,search.median_us/1000,search.worst_us/1000);
        auto cached_search=[&](const char* label,bool empty) {
            cached.clear();
            std::uint64_t evaluations=0;
            const auto result=measure(1,options.rounds,[&] {
                const auto before=cached.source_evaluations();
                const auto found=spp::moon::phases(start,search_end,cached);
                if (found.size()!=events) throw std::runtime_error("cached search returned a different event count");
                double sum=0;
                for (const auto& e:found) sum+=static_cast<double>(e.time.utc().day)+e.time.utc().frac+static_cast<int>(e.phase);
                evaluations=cached.source_evaluations()-before;
                return sum;
            },[&]{if(empty) cached.clear();});
            std::println("  Cached 35-day search ({}), {} events: best {:.3f}, median {:.3f}, worst {:.3f} ms; {} ELP evaluations/search",
                label,events,result.best_us/1000,result.median_us/1000,result.worst_us/1000,evaluations);
        };
        cached_search("repeated workload",false);
        cached_search("empty lunar cache each search",true);
    }
}
#endif

}  // namespace

int benchmark(const Options& options) {
    std::println("compiler {}, CMake config {}, {} build; one thread",compiler_name(),SPP_BENCH_CONFIG,kBuild);
#ifndef NDEBUG
    std::println("Debug build: these timings cannot establish Release performance.");
#endif
#if SPP_BENCH_MOON
    if (options.moon_only) { moon_benchmark(options); return 0; }
#else
    if (options.moon_options)
        throw std::invalid_argument("Moon benchmarks require SPP_BUILD_MOON=ON");
#endif
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
#if SPP_BENCH_MOON
    moon_benchmark(options);
#endif
    return 0;
}

int main(int argc,char** argv) {
    try {
        if (argc==2 && std::string_view(argv[1])=="--help") {
            std::println("usage: spp_bench [--moon-only] [--samples=128] [--rounds=5] [--fit=both|de405|llr]");
            return 0;
        }
        return benchmark(parse_options(argc,argv));
    } catch(const std::exception& e) {
        std::println(stderr,"spp_bench: {}",e.what());
        return 1;
    }
}
