#pragma once

// IAU 1980 nutation theory, truncated to the 20 largest terms
// (Meeus, Astronomical Algorithms, Table 22.A). Truncation error is below
// ~0.05" in Δψ and ~0.02" in Δε.
namespace spp::model::iau1980 {

struct Nutation {
    double dpsi_arcsec;
    double deps_arcsec;
};

// T: Julian centuries of TT from J2000.0.
Nutation nutation(double T);

}  // namespace spp::model::iau1980
