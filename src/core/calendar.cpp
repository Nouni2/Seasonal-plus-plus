#include "spp/core/calendar.hpp"

namespace spp {

bool is_gregorian(const CivilDate& d) {
    if (d.year != 1582) return d.year > 1582;
    if (d.month != 10) return d.month > 10;
    return d.day >= 15;
}

bool is_valid(const CivilDate& d) {
    if (d.month < 1 || d.month > 12 || d.day < 1) return false;
    if (d.year == 1582 && d.month == 10 && d.day >= 5 && d.day <= 14) return false;
    // Floor-safe divisibility for negative (astronomical) years.
    const auto divisible = [](int y, int n) { return y % n == 0; };
    const bool leap = d.year > 1582 ? (divisible(d.year, 4) && (!divisible(d.year, 100) || divisible(d.year, 400)))
                                    : divisible(d.year, 4);
    constexpr int kDays[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    return d.day <= kDays[d.month - 1] + (d.month == 2 && leap ? 1 : 0);
}

}  // namespace spp
