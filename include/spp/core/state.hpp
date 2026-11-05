#pragma once

#include <array>

namespace spp {

// Frame and origin are part of the type: a state from another origin/frame
// requires an explicit transformation before it can be used interchangeably.
namespace origin { struct Geocentric {}; struct Selenocentric {}; struct Relative {}; }
namespace frame {
struct MeanEclipticJ2000 {};
struct Icrf {};
struct MoonMeanEarth {};
}

template<class Origin, class Frame>
struct CartesianPosition {
    std::array<double, 3> position_km;
};

template<class Origin, class Frame>
struct CartesianState {
    std::array<double, 3> position_km;
    std::array<double, 3> velocity_km_per_day;
};

using GeocentricEclipticJ2000State =
    CartesianState<origin::Geocentric, frame::MeanEclipticJ2000>;
using GeocentricEclipticJ2000Position =
    CartesianPosition<origin::Geocentric, frame::MeanEclipticJ2000>;

}  // namespace spp
