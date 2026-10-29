#pragma once

#include <cstdint>
#include <optional>
#include <string>

#include "spp/core/calendar.hpp"
#include "spp/core/constants.hpp"

namespace spp {

// Time-scale tags. An Instant<TT> cannot be passed where an Instant<UT1> is
// expected: mixing scales (a ~69 s error today) is a compile error.
struct UTC {};
struct UT1 {};
struct TT {};

// A Julian Date split into integer day + fraction to keep full precision
// over millennia. frac is in [0, 1) and 0.0 is noon (JD convention).
template <class Scale>
struct Instant {
    std::int64_t day;
    double frac;

    double days_since_j2000() const {
        return static_cast<double>(day - jd_j2000_day) + frac;
    }
    double centuries_since_j2000() const {
        return days_since_j2000() / days_per_century;
    }
};

// An instant carried in UTC, plus the Earth-orientation knowledge (DUT1)
// needed to derive UT1 from it.
//
// DUT1 = UT1 - UTC can be supplied (IERS Bulletin A). When it isn't, it is
// estimated: from 1972 on, where TT - UTC = 32.184 s + (TAI - UTC) is exact,
// as (TT - UTC) - ΔT; before 1972, as 0.
//
// The factories throw std::invalid_argument for a non-finite fraction, second
// or DUT1, and for a civil date that does not exist (is_valid). Hours, minutes
// and seconds outside their usual ranges are accepted and carried into the
// following fields: 24:00 is 0:00 of the next day.
class Time {
public:
    // Normalizes frac into [0, 1), carrying whole days into `day`.
    static Time from_jd(std::int64_t day, double frac, std::optional<double> dut1_s = std::nullopt);
    static Time from_civil(const CivilDateTime& utc, std::optional<double> dut1_s = std::nullopt);
    static Time from_utc(int year, int month, int day,
                         int hour = 0, int minute = 0, double second = 0.0,
                         std::optional<double> dut1_s = std::nullopt);

    Instant<UTC> utc() const;
    Instant<UT1> ut1() const;   // UTC + DUT1
    Instant<TT> tt() const;     // UTC + (TT - UTC); see tt_minus_utc()

    double tt_minus_utc() const;   // seconds: exact from leap seconds when known, else DUT1 + ΔT
    double delta_t() const;        // ΔT = TT - UT1, seconds
    double dut1() const { return dut1_s_; }
    bool dut1_estimated() const { return !dut1_given_; }
    double decimal_year() const;

    Time shifted(double days) const;   // same DUT1 (or re-estimated), renormalized; throws if days is not finite

    CivilDateTime to_civil() const;    // UTC, unrounded seconds
    std::string iso_utc() const;       // iso_date + "THH:MM:SSZ", rounded to the second

private:
    Time(std::int64_t day, double frac, std::optional<double> dut1_s);

    std::int64_t day_;
    double frac_;
    double dut1_s_;
    bool dut1_given_;
};

// "YYYY-MM-DD" (ISO 8601). Years outside 0..9999 use the expanded form with a
// sign and at least four digits: -0100-03-01, +10000-01-01.
std::string iso_date(const CivilDate& d);

// TAI - UTC in seconds at a UTC instant, from the IERS leap-second table, or
// nullopt before 1972-01-01 and after the date up to which the table is known.
std::optional<double> tai_minus_utc(Instant<UTC> t);

}  // namespace spp
