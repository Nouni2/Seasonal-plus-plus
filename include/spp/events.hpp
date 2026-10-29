#pragma once

#include <optional>

#include "spp/core/calendar.hpp"
#include "spp/core/constants.hpp"
#include "spp/core/time.hpp"
#include "spp/core/types.hpp"

// Solar events found by root-finding over time. The functions that solve for
// an event throw std::invalid_argument for a non-finite observer coordinate or
// target altitude.
namespace spp::events {

// Altitude of the Sun's center at each twilight. Sunrise and sunset are not a
// center altitude: they are the instants the refracted upper limb is at 0
// (upper_limb_crossing with target 0, as events::day solves them).
// upper_limb_deg is the center's geometric altitude when the upper limb touches
// the geometric horizon, with the mean semidiameter: use it with refraction off.
inline constexpr double upper_limb_deg = -sun_semidiameter_deg;
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
// `transit` must be the meridian transit for obs.lon_deg, from events::transit
// or Day::transit: the first guess is measured from it, and the altitude there
// is taken as the day's highest. Debug builds assert it (|hour angle| < 0.01°).
std::optional<Time> crossing(const Time& transit, const Observer& obs, const Atmosphere& atm,
                             double target_alt_deg, Direction dir, const SolverConfig& cfg = {});

// Apparent altitude of the Sun's upper limb: the topocentric center raised by
// the semidiameter (959.63" / R at t), then refracted at that altitude.
double upper_limb_altitude_deg(const Time& t, const Observer& obs, const Atmosphere& atm);

// Instant the apparent altitude of the upper limb reaches target_alt_deg on one
// side of transit, or nullopt as for crossing. Target 0 is sunrise or sunset.
// The same precondition on `transit` as crossing.
std::optional<Time> upper_limb_crossing(const Time& transit, const Observer& obs, const Atmosphere& atm,
                                        double target_alt_deg, Direction dir, const SolverConfig& cfg = {});

// Transit, sunrise, sunset and day length for a civil date. Sunrise and sunset
// are the instants the refracted upper limb is at 0; the day type uses the same
// condition at transit and 12 h either side. Next to polar day,
// one crossing can be missing because the Sun is already (or still) up at the
// edge of the 24 h centered on transit; duration_hours is the time above the
// horizon within that window.
Day day(const CivilDate& date, const Observer& obs, const Atmosphere& atm = {},
        std::optional<double> dut1_s = std::nullopt, const SolverConfig& cfg = {});

}  // namespace spp::events
