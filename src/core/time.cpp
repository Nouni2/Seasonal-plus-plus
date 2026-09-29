#include "spp/core/time.hpp"

#include <algorithm>
#include <cmath>
#include <format>

#include "spp/model/espenak_meeus.hpp"
#include "spp/model/meeus.hpp"

namespace spp {

namespace {

// IERS leap-second table (Leap_Second.dat, updated through Bulletin C 72,
// July 2026): TAI - UTC in seconds from the given MJD, 0h UTC, on.
struct LeapStep {
    int mjd;
    int tai_minus_utc;
};

constexpr LeapStep kLeapSeconds[] = {
    {41317, 10}, {41499, 11}, {41683, 12}, {42048, 13}, {42413, 14}, {42778, 15}, {43144, 16},
    {43509, 17}, {43874, 18}, {44239, 19}, {44786, 20}, {45151, 21}, {45516, 22}, {46247, 23},
    {47161, 24}, {47892, 25}, {48257, 26}, {48804, 27}, {49169, 28}, {49534, 29}, {50083, 30},
    {50630, 31}, {51179, 32}, {53736, 33}, {54832, 34}, {56109, 35}, {57204, 36}, {57754, 37},
};

// The table is known to hold until this date (MJD 61584 = 2027-06-28, 0h UTC).
constexpr std::int64_t kLeapSecondsValidUntilMjd = 61584;

// TT - TAI, seconds.
constexpr double kTtMinusTai = 32.184;

// Integer MJD of the UTC civil day containing the instant (day, frac).
std::int64_t mjd_day(std::int64_t day, double frac) { return frac >= 0.5 ? day - 2400000 : day - 2400001; }

// Brings frac into [0, 1), carrying whole days into day.
void normalize(std::int64_t& day, double& frac) {
    const double whole = std::floor(frac);
    day += static_cast<std::int64_t>(whole);
    frac -= whole;
    if (frac >= 1.0) {   // e.g. -1e-17 + 1.0 rounds to 1.0
        day += 1;
        frac -= 1.0;
    }
}

template <class Scale>
Instant<Scale> offset(std::int64_t day, double frac, double seconds) {
    frac += seconds / seconds_per_day;
    normalize(day, frac);
    return {day, frac};
}

double decimal_year_of(std::int64_t day, double frac) {
    return 2000.0 + Instant<UTC>{day, frac}.centuries_since_j2000() * 100.0;
}

double delta_t_model(std::int64_t day, double frac) {
    return model::espenak_meeus::delta_t_seconds(decimal_year_of(day, frac));
}

std::optional<double> tai_minus_utc_at(std::int64_t day, double frac) {
    const std::int64_t mjd = mjd_day(day, frac);
    if (mjd < kLeapSeconds[0].mjd || mjd >= kLeapSecondsValidUntilMjd) return std::nullopt;
    const auto* step = std::upper_bound(std::begin(kLeapSeconds), std::end(kLeapSeconds), mjd,
                                        [](std::int64_t m, const LeapStep& s) { return m < s.mjd; });
    return static_cast<double>((step - 1)->tai_minus_utc);
}

// DUT1 when none is supplied. Where TT - UTC is exact, UT1 - UTC follows from
// the ΔT model; it is clamped to the ±0.9 s that the leap-second system
// maintains. After the table, the last value is held, which keeps TT
// continuous; before 1972, UTC is taken as UT1.
double estimated_dut1(std::int64_t day, double frac) {
    const std::int64_t mjd = mjd_day(day, frac);
    if (mjd < kLeapSeconds[0].mjd) return 0.0;
    if (mjd >= kLeapSecondsValidUntilMjd) {
        day = kLeapSecondsValidUntilMjd + 2400000 - 1;   // noon of the last day covered
        frac = 0.0;
    }
    const double tt_minus_utc = kTtMinusTai + *tai_minus_utc_at(day, frac);
    return std::clamp(tt_minus_utc - delta_t_model(day, frac), -0.9, 0.9);
}

// Floor division derived from fmod, consistent with floor_mod below at exact
// multiples (std::floor(x / y) can be off by one there after rounding).
double floor_div(double x, double y) {
    double mod = std::fmod(x, y);
    double div = (x - mod) / y;
    if (mod != 0.0 && ((y < 0.0) != (mod < 0.0))) div -= 1.0;
    if (div == 0.0) return std::copysign(0.0, x / y);
    double floordiv = std::floor(div);
    if (div - floordiv > 0.5) floordiv += 1.0;
    return floordiv;
}

double floor_mod(double x, double y) {
    double mod = std::fmod(x, y);
    if (mod != 0.0) {
        if ((y < 0.0) != (mod < 0.0)) mod += y;
    } else {
        mod = std::copysign(0.0, y);
    }
    return mod;
}

}  // namespace

std::optional<double> tai_minus_utc(Instant<UTC> t) { return tai_minus_utc_at(t.day, t.frac); }

Time::Time(std::int64_t day, double frac, std::optional<double> dut1_s)
    : day_(day), frac_(frac), dut1_s_(0.0), dut1_given_(dut1_s.has_value()) {
    normalize(day_, frac_);
    dut1_s_ = dut1_s ? *dut1_s : estimated_dut1(day_, frac_);
}

Time Time::from_jd(std::int64_t day, double frac, std::optional<double> dut1_s) { return Time(day, frac, dut1_s); }

Time Time::from_civil(const CivilDateTime& utc, std::optional<double> dut1_s) {
    const auto jd = model::meeus::jd_from_civil(utc);
    return Time(jd.day, jd.frac, dut1_s);
}

Time Time::from_utc(int year, int month, int day, int hour, int minute, double second,
                    std::optional<double> dut1_s) {
    return from_civil({year, month, day, hour, minute, second}, dut1_s);
}

Instant<UTC> Time::utc() const { return {day_, frac_}; }

Instant<UT1> Time::ut1() const {
    if (dut1_s_ == 0.0) return {day_, frac_};
    return offset<UT1>(day_, frac_, dut1_s_);
}

double Time::tt_minus_utc() const {
    if (const auto tai = tai_minus_utc_at(day_, frac_)) return kTtMinusTai + *tai;
    return dut1_s_ + delta_t_model(day_, frac_);
}

Instant<TT> Time::tt() const { return offset<TT>(day_, frac_, tt_minus_utc()); }

double Time::decimal_year() const { return decimal_year_of(day_, frac_); }

double Time::delta_t() const { return tt_minus_utc() - dut1_s_; }

Time Time::shifted(double days) const {
    // Whole days go to the integer part, so large shifts keep full precision.
    const double whole = std::floor(days);
    return Time(day_ + static_cast<std::int64_t>(whole), frac_ + (days - whole),
                dut1_given_ ? std::optional<double>(dut1_s_) : std::nullopt);
}

CivilDateTime Time::to_civil() const {
    const auto c = model::meeus::calendar_from_jd(day_, frac_);
    const double total_seconds = c.day_fraction * seconds_per_day;
    const int hour = static_cast<int>(floor_div(total_seconds, 3600.0));
    const int minute = static_cast<int>(floor_div(floor_mod(total_seconds, 3600.0), 60.0));
    const double second = total_seconds - hour * 3600 - minute * 60;
    return {c.year, c.month, c.day, hour, minute, second};
}

std::string Time::iso_utc() const {
    // Round at the JD level so 23:59:59.7 correctly rolls over the date.
    std::int64_t day = day_;
    double frac = frac_ + 0.5 / seconds_per_day;
    normalize(day, frac);
    const auto c = model::meeus::calendar_from_jd(day, frac);
    const auto total = static_cast<int>(c.day_fraction * seconds_per_day);
    return std::format("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}Z",
                       c.year, c.month, c.day, total / 3600, (total % 3600) / 60, total % 60);
}

}  // namespace spp
