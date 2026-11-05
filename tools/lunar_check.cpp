#include "spp/moon.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "error_stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <format>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

using spp::model::elp_mpp02::Ephemeris;
using spp::model::elp_mpp02::Fit;
static_assert(!std::is_invocable_v<Ephemeris, spp::Instant<spp::TT>>);
static_assert(spp::moon::Ephemeris<const Ephemeris>);
static_assert(spp::moon::PositionEphemeris<const Ephemeris>);
static_assert(!std::is_invocable_v<decltype(&Ephemeris::position), const Ephemeris&, spp::Instant<spp::TT>>);

// A stateful alternative demonstrates body/model separation without changing
// the Moon API or requiring all providers to be immutable analytical models.
struct OtherModel {
    int calls = 0;
    spp::moon::State operator()(spp::Instant<spp::TDB>) {
        ++calls;
        return {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    }
};

struct WrongFrameModel {
    spp::CartesianState<spp::origin::Geocentric, void> operator()(spp::Instant<spp::TDB>);
};
static_assert(!spp::moon::Ephemeris<WrongFrameModel>);

struct WrongPositionFrame {
    spp::CartesianPosition<spp::origin::Geocentric, void> position(spp::Instant<spp::TDB>);
};
static_assert(!spp::moon::PositionSource<WrongPositionFrame>);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Exception, class F>
void throws(F function) {
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("expected exception was not raised");
}

void contract_checks() {
    const Ephemeris model;
    const spp::Instant<spp::TDB> epoch{2451545, 0.25};
    const auto baseline = model(epoch);
    for (const auto split : {spp::Instant<spp::TDB>{2451544, 1.25},
                             spp::Instant<spp::TDB>{2451546, -0.75}}) {
        const auto result = model(split);
        require(result.position_km == baseline.position_km
                && result.velocity_km_per_day == baseline.velocity_km_per_day,
                "equivalent split dates differ");
    }
    const auto default_state = spp::moon::geometric(epoch);
    require(default_state.position_km == baseline.position_km
            && default_state.velocity_km_per_day == baseline.velocity_km_per_day,
            "default Moon model is not full ELP/MPP02 DE405");
    require(spp::moon::position(epoch).position_km == model.position(epoch).position_km,
            "default Moon position is not full ELP/MPP02 DE405");
    OtherModel other;
    const auto alternate = spp::moon::geometric(epoch, other);
    require(other.calls == 1 && alternate.position_km[0] == 1.0,
            "explicit stateful model was not used");
    const auto alternate_position = spp::moon::position(epoch, other);
    require(other.calls == 2 && alternate_position.position_km == alternate.position_km,
            "full-state fallback for position was not used");
    const Ephemeris llr(Fit::LLR);
    require(llr(epoch).position_km != baseline.position_km, "fitted sets are indistinguishable");
    throws<std::invalid_argument>([] { Ephemeris invalid(static_cast<Fit>(42)); });
    for (const double fraction : {std::numeric_limits<double>::infinity(),
                                  -std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::quiet_NaN()})
    {
        throws<std::invalid_argument>([&] { model({2451545, fraction}); });
        throws<std::invalid_argument>([&] { model.position({2451545, fraction}); });
    }
    for (const auto invalid : {spp::Instant<spp::TDB>{std::numeric_limits<std::int64_t>::max(), 0.0},
                               spp::Instant<spp::TDB>{std::numeric_limits<std::int64_t>::min(), 0.0},
                               spp::Instant<spp::TDB>{2451545, 1e100},
                               spp::Instant<spp::TDB>{1319270, -0.25},
                               spp::Instant<spp::TDB>{2853320, 0.25}})
    {
        throws<std::out_of_range>([&] { model(invalid); });
        throws<std::out_of_range>([&] { model.position(invalid); });
    }
    for (const auto split : {spp::Instant<spp::TDB>{2451544, 1.25},
                             spp::Instant<spp::TDB>{2451546, -0.75}})
        require(model.position(split).position_km == model.position(epoch).position_km,
                "equivalent split position dates differ");

    // Five-point difference reduces truncation error while allowing a larger
    // step: at ancient epochs, tiny steps amplify double argument-rounding
    // noise. This checks the complete position derivative independently.
    constexpr double step = 1.0 / 32.0;
    double max_difference = 0.0;
    for (const Fit fit : {Fit::DE405, Fit::LLR}) {
        const Ephemeris provider(fit);
        for (const std::int64_t day : {1356000LL, 2451545LL, 2816000LL}) {
            const auto center = provider({day, 0.5});
            const auto before = provider({day, 0.5 - step});
            const auto after = provider({day, 0.5 + step});
            const auto before2 = provider({day, 0.5 - 2*step});
            const auto after2 = provider({day, 0.5 + 2*step});
            for (std::size_t i = 0; i < 3; ++i)
                max_difference = std::max(max_difference,
                    std::abs((before2.position_km[i] - 8*before.position_km[i]
                              + 8*after.position_km[i] - after2.position_km[i]) / (12*step)
                             - center.velocity_km_per_day[i]));
        }
    }
    std::printf("Five-point velocity difference: %.9g km/day\n", max_difference);
    require(max_difference < 0.002, "analytical velocity fails finite-difference check");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: spp_lunar_check <reference.csv>");
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("cannot open independent reference states");
        std::string line;
        int count = 0;
        double max_position = 0.0, max_velocity = 0.0;
        spp::check::ErrorStats states, positions, velocities, agreement;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream row(line);
            int fit = -1;
            spp::Instant<spp::TDB> time{};
            spp::moon::State expected{};
            require(static_cast<bool>(row >> fit >> time.day >> time.frac), "invalid reference epoch");
            require(fit == 0 || fit == 1, "invalid reference fit");
            for (double& value : expected.position_km)
                require(static_cast<bool>(row >> value) && std::isfinite(value), "invalid reference position");
            for (double& value : expected.velocity_km_per_day)
                require(static_cast<bool>(row >> value) && std::isfinite(value), "invalid reference velocity");
            std::string extra;
            require(!(row >> extra), "extra reference column");
            const Ephemeris model(fit == 0 ? Fit::LLR : Fit::DE405);
            const auto state = spp::moon::geometric(time, model);
            const auto position = spp::moon::position(time, model);
            const auto context = std::format("{} TDB JD {}+{:.12f}",
                fit == 0 ? "LLR" : "DE405", time.day, time.frac);
            for (std::size_t i = 0; i < 3; ++i) {
                require(std::isfinite(state.position_km[i]) && std::isfinite(state.velocity_km_per_day[i]),
                        "non-finite lunar state");
                max_position = std::max(max_position, std::abs(state.position_km[i] - expected.position_km[i]));
                max_position = std::max(max_position, std::abs(position.position_km[i] - expected.position_km[i]));
                // Both branches evaluate every coefficient. Allow libm's last
                // bits to differ when the compiler selects sincos for the
                // state and a single sin/cos for position-only evaluation.
                require(std::isfinite(position.position_km[i]), "non-finite position-only result");
                max_velocity = std::max(max_velocity, std::abs(state.velocity_km_per_day[i] - expected.velocity_km_per_day[i]));
                const auto component = std::format("{} component {}", context, "xyz"[i]);
                states.add((state.position_km[i]-expected.position_km[i])*1000.0, component);
                positions.add((position.position_km[i]-expected.position_km[i])*1000.0, component);
                velocities.add(state.velocity_km_per_day[i]-expected.velocity_km_per_day[i], component);
                agreement.add((position.position_km[i]-state.position_km[i])*1000.0, component);
            }
            ++count;
        }
        require(!input.bad(), "reference read failure");
        require(count >= 100, "reference coverage is incomplete");
        std::printf("ELP/MPP02: %d states, max component error %.9g km, %.9g km/day\n",
                    count, max_position, max_velocity);
        std::puts("Implementation agreement with pinned upstream ELP/MPP02 (both fits):");
        const bool states_ok = states.report("full-state position", "m", 1e-4);
        const bool positions_ok = positions.report("position-only", "m", 1e-4);
        const bool velocities_ok = velocities.report("analytic velocity", "km/day", 1e-7);
        const bool agreement_ok = agreement.report("position vs full state", "m", 1e-4);
        require(states_ok && positions_ok && velocities_ok && agreement_ok, "ELP component residual limits exceeded");
        std::puts("These references check implementation agreement, not independent physical accuracy.");
        require(max_position < 1e-7 && max_velocity < 1e-7, "reference comparison failed");
        contract_checks();
        std::puts("PASS: references, velocity, split dates, invalid inputs and model selection");
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "FAIL: %s\n", error.what());
        return 1;
    }
}
