#include "spp/core/time.hpp"

#include <cmath>
#include <format>

#include "spp/model/espenak_meeus.hpp"
#include "spp/model/meeus.hpp"

namespace spp {

namespace {

// Carries whole days out of a fraction that left [0, 1) after adding an offset.
template <class Scale>
Instant<Scale> renormalized(std::int64_t day, double frac) {
    if (frac >= 1.0) {
        const auto extra = static_cast<std::int64_t>(frac);
        day += extra;
        frac -= static_cast<double>(extra);
    } else if (frac < 0.0) {
        const auto extra = static_cast<std::int64_t>(std::ceil(-frac));
        day -= extra;
        frac += static_cast<double>(extra);
    }
    return {day, frac};
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

Time::Time(std::int64_t day, double frac, double dut1_s) : dut1_s_(dut1_s) {
    // Truncate toward zero, then pull a negative fraction back into [0, 1).
    const auto extra = static_cast<std::int64_t>(frac);
    day_ = day + extra;
    frac_ = frac - static_cast<double>(extra);
    if (frac_ < 0.0) {
        day_ -= 1;
        frac_ += 1.0;
    }
}

Time Time::from_jd(std::int64_t day, double frac, double dut1_s) {
    return Time(day, frac, dut1_s);
}

Time Time::from_civil(const CivilDateTime& utc, double dut1_s) {
    const auto jd = model::meeus::jd_from_civil(utc);
    return Time(jd.day, jd.frac, dut1_s);
}

Time Time::from_utc(int year, int month, int day, int hour, int minute, double second, double dut1_s) {
    return from_civil({year, month, day, hour, minute, second}, dut1_s);
}

Instant<UTC> Time::utc() const { return {day_, frac_}; }

Instant<UT1> Time::ut1() const {
    if (dut1_s_ == 0.0) return {day_, frac_};
    return renormalized<UT1>(day_, frac_ + dut1_s_ / seconds_per_day);
}

Instant<TT> Time::tt() const {
    const double dt_seconds = delta_t() + dut1_s_;
    return renormalized<TT>(day_, frac_ + dt_seconds / seconds_per_day);
}

double Time::decimal_year() const {
    // Approximate; only used to select and evaluate the Delta T fits.
    return 2000.0 + utc().centuries_since_j2000() * 100.0;
}

double Time::delta_t() const {
    return model::espenak_meeus::delta_t_seconds(decimal_year());
}

Time Time::shifted(double days) const { return Time(day_, frac_ + days, dut1_s_); }

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
    if (frac >= 1.0) {
        day += 1;
        frac -= 1.0;
    }
    const auto c = model::meeus::calendar_from_jd(day, frac);
    const auto total = static_cast<int>(c.day_fraction * seconds_per_day);
    return std::format("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}Z",
                       c.year, c.month, c.day, total / 3600, (total % 3600) / 60, total % 60);
}

}  // namespace spp
