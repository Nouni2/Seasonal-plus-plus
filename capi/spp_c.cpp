// The C interface (include/spp/spp_c.h): converts between C types and the
// library's, validates inputs, and turns every C++ exception into a status
// code. It calls the public API only and computes no astronomy of its own.

#include "spp/spp_c.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>

#include "spp/spp.hpp"

#ifndef SPP_C_VERSION_STRING
#define SPP_C_VERSION_STRING "unknown"
#endif

// The layout the header documents, checked where the library is compiled.
static_assert(sizeof(spp_observer) == 24, "spp_observer layout");
static_assert(sizeof(spp_solar_day) == 72, "spp_solar_day layout");
static_assert(offsetof(spp_solar_day, transit_utc_ms) == 24, "spp_solar_day layout");
static_assert(offsetof(spp_solar_day, duration_hours) == 64, "spp_solar_day layout");

namespace {

constexpr std::int64_t kMsPerDay = 86'400'000;
constexpr std::int64_t kUnixEpochJdDay = 2440587;   // 1970-01-01T00:00Z is JD 2440587.5
constexpr int kFirstYear = -1000;
constexpr int kLastYear = 3000;
constexpr std::int32_t kMaxCount = 1'000'000;

std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    const std::int64_t q = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? q - 1 : q;
}

// UTC milliseconds -> Time (DUT1 estimated): the JD split is built from the
// whole days and the millisecond remainder, so no large value is rounded.
spp::Time from_utc_ms(std::int64_t ms) {
    const std::int64_t days = floor_div(ms, kMsPerDay);
    const std::int64_t rest = ms - days * kMsPerDay;   // [0, 86 400 000)
    return spp::Time::from_jd(kUnixEpochJdDay + days, 0.5 + static_cast<double>(rest) / kMsPerDay);
}

// Time -> UTC milliseconds, rounded to the nearest millisecond.
std::int64_t to_utc_ms(const spp::Time& t) {
    const auto u = t.utc();
    return (u.day - kUnixEpochJdDay) * kMsPerDay + std::llround((u.frac - 0.5) * static_cast<double>(kMsPerDay));
}

// The supported span, [-1000-01-01T00:00Z, 3001-01-01T00:00Z).
std::int64_t first_ms() {
    static const std::int64_t v = to_utc_ms(spp::Time::from_utc(kFirstYear, 1, 1));
    return v;
}
std::int64_t end_ms() {
    static const std::int64_t v = to_utc_ms(spp::Time::from_utc(kLastYear + 1, 1, 1));
    return v;
}
bool in_range(std::int64_t ms) { return ms >= first_ms() && ms < end_ms(); }

bool valid(const spp_observer* o) {
    return o != nullptr && std::isfinite(o->lat_deg) && std::isfinite(o->lon_deg) && std::isfinite(o->elev_m) &&
           o->lat_deg >= -90.0 && o->lat_deg <= 90.0 && o->lon_deg >= -180.0 && o->lon_deg <= 180.0 &&
           o->elev_m >= -1000.0 && o->elev_m <= 11000.0;
}

spp::Observer observer(const spp_observer* o) { return {o->lat_deg, o->lon_deg, o->elev_m}; }

// Runs `body`, mapping exceptions to status codes.
template <class F>
std::int32_t guarded(F body) {
    try {
        return body();
    } catch (const std::invalid_argument&) {
        return SPP_E_ARGUMENT;
    } catch (...) {
        return SPP_E_INTERNAL;
    }
}

void put_event(const std::optional<spp::Time>& t, std::int32_t& has, std::int64_t& ms) {
    has = t ? 1 : 0;
    ms = t ? to_utc_ms(*t) : 0;
}

}  // namespace

extern "C" {

SPP_C_API std::int32_t SPP_C_CALL spp_c_abi_version(void) { return SPP_C_ABI_VERSION; }

SPP_C_API const char* SPP_C_CALL spp_c_library_version(void) { return SPP_C_VERSION_STRING; }

SPP_C_API std::int32_t SPP_C_CALL spp_c_solar_day_size(void) { return static_cast<std::int32_t>(sizeof(spp_solar_day)); }

SPP_C_API std::int32_t SPP_C_CALL spp_solar_day_utc(const spp_observer* obs, std::int32_t year, std::int32_t month,
                                                     std::int32_t day, spp_solar_day* out) {
    if (out == nullptr || out->struct_size != static_cast<std::int32_t>(sizeof(spp_solar_day)) || !valid(obs))
        return SPP_E_ARGUMENT;
    if (year < kFirstYear || year > kLastYear) return SPP_E_RANGE;
    if (!spp::is_valid({year, month, day})) return SPP_E_ARGUMENT;
    return guarded([&] {
        const spp::Observer o = observer(obs);
        const auto d = spp::events::day({year, month, day}, o);
        using spp::events::Direction;
        const auto dawn = spp::events::crossing(d.transit, o, {}, spp::events::civil_twilight_deg, Direction::Before);
        const auto dusk = spp::events::crossing(d.transit, o, {}, spp::events::civil_twilight_deg, Direction::After);

        spp_solar_day r{};
        r.struct_size = out->struct_size;
        r.day_type = d.type == spp::events::DayType::PolarDay     ? SPP_DAY_POLAR_DAY
                     : d.type == spp::events::DayType::PolarNight ? SPP_DAY_POLAR_NIGHT
                                                                  : SPP_DAY_NORMAL;
        r.transit_utc_ms = to_utc_ms(d.transit);
        put_event(d.rise, r.has_rise, r.rise_utc_ms);
        put_event(d.set, r.has_set, r.set_utc_ms);
        put_event(dawn, r.has_civil_dawn, r.civil_dawn_utc_ms);
        put_event(dusk, r.has_civil_dusk, r.civil_dusk_utc_ms);
        r.duration_hours = d.duration_hours;
        *out = r;
        return SPP_OK;
    });
}

SPP_C_API std::int32_t SPP_C_CALL spp_upper_limb_altitude(const spp_observer* obs, std::int64_t utc_ms,
                                                          double* out_alt_deg) {
    if (out_alt_deg == nullptr || !valid(obs)) return SPP_E_ARGUMENT;
    if (!in_range(utc_ms)) return SPP_E_RANGE;
    return guarded([&] {
        *out_alt_deg = spp::events::upper_limb_altitude_deg(from_utc_ms(utc_ms), observer(obs), {});
        return SPP_OK;
    });
}

SPP_C_API std::int32_t SPP_C_CALL spp_upper_limb_altitudes(const spp_observer* obs, std::int64_t start_utc_ms,
                                                           std::int64_t step_ms, std::int32_t count,
                                                           double* out_alt_deg) {
    if (!valid(obs) || count < 0 || count > kMaxCount || (count > 0 && out_alt_deg == nullptr))
        return SPP_E_ARGUMENT;
    if (count == 0) return SPP_OK;
    if (!in_range(start_utc_ms)) return SPP_E_RANGE;
    // The last instant must be in range; bound the step first so the product cannot overflow.
    const std::int64_t span = end_ms() - first_ms();
    if (count > 1 && (step_ms > span / (count - 1) || step_ms < -span / (count - 1))) return SPP_E_RANGE;
    if (!in_range(start_utc_ms + static_cast<std::int64_t>(count - 1) * step_ms)) return SPP_E_RANGE;
    return guarded([&] {
        const spp::Observer o = observer(obs);
        for (std::int32_t i = 0; i < count; ++i)
            out_alt_deg[i] =
                spp::events::upper_limb_altitude_deg(from_utc_ms(start_utc_ms + i * step_ms), o, {});
        return SPP_OK;
    });
}

}  // extern "C"
