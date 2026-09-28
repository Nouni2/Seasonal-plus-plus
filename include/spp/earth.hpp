#pragma once

#include "spp/core/time.hpp"

// Earth orientation: rotation (sidereal time), nutation, obliquity, ΔT.
namespace spp::earth {

struct Nutation {
    double dpsi_deg;             // Δψ, nutation in longitude
    double deps_deg;             // Δε, nutation in obliquity
    double mean_obliquity_deg;   // ε0
    double true_obliquity_deg;   // ε = ε0 + Δε
};

Nutation nutation(Instant<TT> t);
double mean_obliquity_deg(Instant<TT> t);

double delta_t(const Time& t);   // TT - UT1, seconds

double gmst_deg(Instant<UT1> t);                        // Greenwich mean sidereal time
double equation_of_equinoxes_deg(const Nutation& n);   // Δψ cos ε
double gast_deg(Instant<UT1> t, const Nutation& n);    // GMST + equation of the equinoxes
double gast_deg(const Time& t);                         // computes nutation at t.tt()

// Local apparent sidereal time, degrees [0, 360). Longitude positive East.
double lst_deg(Instant<UT1> t, const Nutation& n, double lon_east_deg);
double lst_deg(const Time& t, double lon_east_deg);

}  // namespace spp::earth
