#pragma once

// Laskar (1986) mean obliquity of the ecliptic (Meeus 22.3),
// valid for about ±10,000 years around J2000.0.
namespace spp::model::laskar {

// T: Julian centuries of TT from J2000.0.
double mean_obliquity_deg(double T);

}  // namespace spp::model::laskar
