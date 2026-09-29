#include "spp/model/meeus.hpp"

#include <cmath>

#include "spp/core/angles.hpp"
#include "spp/core/constants.hpp"

namespace spp::model::meeus {

SplitJd jd_from_civil(const CivilDateTime& c) {
    // Algorithm 7.1: January/February count as months 13/14 of the previous year.
    int year = c.year;
    int month = c.month;
    if (month <= 2) {
        year -= 1;
        month += 12;
    }

    // The Gregorian correction B applies from 1582-10-15; Julian calendar before.
    double B = 0.0;
    if (is_gregorian(c.date())) {
        const double A = std::floor(year / 100.0);
        B = 2.0 - A + std::floor(A / 4.0);
    }

    const double jd_noon = std::floor(365.25 * (year + 4716)) + std::floor(30.6001 * (month + 1))
                           + c.day + B - 1524.0;

    // JD starts at noon: civil midnight is JD .5
    const double day_fraction = (c.hour / 24.0) + (c.minute / 1440.0) + (c.second / seconds_per_day);
    return {static_cast<std::int64_t>(jd_noon), day_fraction - 0.5};
}

CalendarDay calendar_from_jd(std::int64_t day, double frac) {
    // Shift by half a day so the split lines up with civil midnight.
    const double g = frac + 0.5;
    const std::int64_t Z = day + static_cast<std::int64_t>(std::floor(g));
    const double F = g - std::floor(g);

    // Floors throughout (not truncation), so negative Julian Days also convert.
    const auto ifloor = [](double x) { return static_cast<std::int64_t>(std::floor(x)); };

    std::int64_t A = Z;
    if (Z >= jd_gregorian_start) {
        const std::int64_t alpha = ifloor((static_cast<double>(Z) - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - ifloor(static_cast<double>(alpha) / 4.0);
    }

    const std::int64_t B = A + 1524;
    const std::int64_t C = ifloor((static_cast<double>(B) - 122.1) / 365.25);
    const std::int64_t D = ifloor(365.25 * static_cast<double>(C));
    const std::int64_t E = ifloor(static_cast<double>(B - D) / 30.6001);

    const double day_with_frac = static_cast<double>(B - D - ifloor(30.6001 * static_cast<double>(E))) + F;
    const double day_floor = std::floor(day_with_frac);

    const int month = static_cast<int>(E < 14 ? E - 1 : E - 13);
    const int year = static_cast<int>(month > 2 ? C - 4716 : C - 4715);
    return {year, month, static_cast<int>(day_floor), day_with_frac - day_floor};
}

double gmst_deg(std::int64_t days, double frac) {
    // Formula 12.4 (IAU 1982); D = days + frac and T are counted on UT1:
    // θ0 = 280.46061837 + 360.98564736629 D + 0.000387933 T² - T³/38710000.
    constexpr double a = 280.46061837;
    constexpr double c = 0.000387933;
    constexpr double d = -1.0 / 38710000.0;

    // 360.98564736629 D = 360 D + r D. With D = days + frac, 360·days vanishes
    // mod 360. r is split into a 24-bit head (so head·days is exact for any
    // |days| < 2^29) and a small tail.
    constexpr double r = 0.98564736629;
    constexpr double r_head = static_cast<double>(static_cast<float>(r));
    constexpr double r_tail = r - r_head;

    const double day_d = static_cast<double>(days);
    const double whole = std::fmod(r_head * day_d, 360.0);   // exact
    const double T = (day_d + frac) / days_per_century;
    return norm360(a + 360.0 * frac + whole + r_tail * day_d + r * frac + c * T * T + d * T * T * T);
}

}  // namespace spp::model::meeus
