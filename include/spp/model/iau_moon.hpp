#pragma once

#include "spp/core/rotation.hpp"
#include "spp/core/time.hpp"

namespace spp::model::iau_moon {

struct Orientation {
    double pole_ra_deg,pole_dec_deg,prime_meridian_deg;
    double pole_ra_rate_deg_day,pole_dec_rate_deg_day,prime_meridian_rate_deg_day;
    FrameRotation<frame::Icrf,frame::MoonMeanEarth> rotation;
};

// IAU2009 trigonometric approximation of the lunar Mean Earth/Polar Axis
// frame, as published in NAIF pck00011.tpc (IAU_MOON). This is a separate
// rotation provider; it does not evaluate or select an orbital ephemeris.
// High-precision DE-series lunar librations require a different provider.
class Rotation {
public:
    static constexpr Instant<TDB> first_epoch{1319270,0};
    static constexpr Instant<TDB> last_epoch{2853320,0};
    // These numerical bounds match the current lunar engine; they are not
    // a claim of physical orientation accuracy over that entire interval.
    Orientation operator()(Instant<TDB> t) const;
};
}  // namespace spp::model::iau_moon
