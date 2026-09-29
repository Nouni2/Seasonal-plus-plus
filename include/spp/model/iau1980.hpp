#pragma once

// IAU 1980 nutation theory, full series (106 terms; Seidelmann 1992,
// Explanatory Supplement to the Astronomical Almanac).
namespace spp::model::iau1980 {

struct Nutation {
    double dpsi_arcsec;
    double deps_arcsec;
};

// T: Julian centuries of TT from J2000.0.
Nutation nutation(double T);

}  // namespace spp::model::iau1980
