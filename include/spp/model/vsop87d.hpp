#pragma once

// VSOP87 version D (Bretagnon & Francou 1988): heliocentric ecliptic
// coordinates referred to the dynamical ecliptic and equinox of date.
namespace spp::model::vsop87d {

struct Term {
    double a, b, c;   // a * cos(b + c * t)
};

struct Spherical {
    double l_rad;   // longitude, not normalized
    double b_rad;   // latitude
    double r_au;    // radius vector
};

// t: Julian millennia of TT from J2000.0, i.e. T / 10.
Spherical earth(double t);

}  // namespace spp::model::vsop87d
