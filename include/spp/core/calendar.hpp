#pragma once

namespace spp {

// Civil calendar date. Julian calendar before 1582-10-15, Gregorian on/after.
// Astronomical year numbering: 1 BC is year 0, 2 BC is year -1.
struct CivilDate {
    int year;
    int month;   // 1..12
    int day;     // 1..31
};

struct CivilDateTime {
    int year;
    int month;
    int day;
    int hour = 0;
    int minute = 0;
    double second = 0.0;

    CivilDate date() const { return {year, month, day}; }
};

bool is_gregorian(const CivilDate& d);   // on/after 1582-10-15

}  // namespace spp
