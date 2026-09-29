#pragma once

#include <optional>

#include "spp/core/calendar.hpp"
#include "spp/core/constants.hpp"
#include "spp/core/time.hpp"
#include "spp/core/types.hpp"

// Solar events found by root-finding over time.
namespace spp::events {

// Apparent altitude of the Sun's center at each event. events::day uses the
// actual semidiameter (sun::semidiameter_deg); upper_limb_deg is its mean.
inline constexpr double upper_limb_deg = -sun_semidiameter_deg;   // sunrise / sunset
inline constexpr double civil_twilight_deg = -6.0;
inline constexpr double nautical_twilight_deg = -12.0;
inline constexpr double astronomical_twilight_deg = -18.0;

enum class Direction { Before = -1, After = 1 };   // rise side / set side of transit

enum class DayType { Normal, PolarDay, PolarNight };

// Convergence of transit and crossing refinement.
struct SolverConfig {
    double time_tolerance_s = 1e-6;
    int max_iter = 60;
};

struct Day {
    CivilDate date;
    DayType type;
    Time transit;
    std::optional<Time> rise;
    std::optional<Time> set;
    double duration_hours;
};

// Full pipeline: apparent, topocentric, refracted altitude of the Sun's center.
double altitude_deg(const Time& t, const Observer& obs, const Atmosphere& atm);

// Local meridian transit (hour angle of the apparent geocentric Sun = 0),
// refined from an approximate instant.
Time transit(const Time& approx, double lon_east_deg, const SolverConfig& cfg = {});

// Instant the Sun's center reaches target_alt_deg on one side of transit,
// or nullopt if it never does within 12 h (or refinement fails).
std::optional<Time> crossing(const Time& transit, const Observer& obs, const Atmosphere& atm,
                             double target_alt_deg, Direction dir, const SolverConfig& cfg = {});

// Transit, sunrise, sunset and day length for a civil date. Next to polar day,
// one crossing can be missing because the Sun is already (or still) up at the
// edge of the 24 h centered on transit; duration_hours is the time above the
// horizon within that window.
Day day(const CivilDate& date, const Observer& obs, const Atmosphere& atm = {},
        std::optional<double> dut1_s = std::nullopt, const SolverConfig& cfg = {});

}  // namespace spp::events
