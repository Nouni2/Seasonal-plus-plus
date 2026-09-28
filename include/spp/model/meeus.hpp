#pragma once

#include <cstdint>

#include "spp/core/calendar.hpp"

// Algorithms from Meeus, Astronomical Algorithms (2nd ed.).
namespace spp::model::meeus {

struct SplitJd {
    std::int64_t day;   // integer JD
    double frac;        // not normalized; 0.0 is noon
};

struct CalendarDay {
    int year;
    int month;
    int day;
    double day_fraction;   // time of day from civil midnight, [0, 1)
};

// Chapter 7: civil date <-> Julian Date, Julian/Gregorian switch at 1582-10-15.
SplitJd jd_from_civil(const CivilDateTime& c);
CalendarDay calendar_from_jd(std::int64_t day, double frac);

// Formula 12.4: GMST in degrees [0, 360). Days/centuries counted on UT1.
double gmst_deg(double days_since_j2000_ut1);

// Same formula evaluated without the ~1e8 deg intermediate of 360.98...·D
// (which a double resolves only to ~1e-8 deg): D is split into whole days and
// a fraction, and 360·(whole days) ≡ 0 (mod 360) is dropped exactly.
double gmst_deg(std::int64_t whole_days_since_j2000_ut1, double frac);

}  // namespace spp::model::meeus
