#pragma once

#include "vector.hpp"
#include "spp/earth.hpp"
#include "spp/core/rotation.hpp"

namespace spp::detail {

// Dynamical mean ecliptic/equinox J2000 <-> equatorial of date. These are
// rotations, not ICRF frame-bias adapters. Precession: IAU 2006 P03;
// nutation and obliquity: Spp's existing IAU 1980 / Laskar conventions.
Vector ecliptic_j2000_to_equatorial(Vector v, Instant<TT> t, const earth::Nutation& n);
Vector equatorial_to_ecliptic_j2000(Vector v, Instant<TT> t, const earth::Nutation& n);
Vector ecliptic_j2000_to_mean_equatorial(Vector v, Instant<TT> t);
Vector mean_equatorial_to_ecliptic_j2000(Vector v, Instant<TT> t);

// Basis-only rotations for repeated vectors at one observation epoch.
// Physical vectors use these directly; coordinate-state derivatives belong
// to FrameRotation, not to this instantaneous observer basis.
struct EquatorialFrame {
    RotationMatrix forward;
    Vector to_equatorial(Vector v) const {
        return {dot(forward[0],v),dot(forward[1],v),dot(forward[2],v)};
    }
    Vector to_ecliptic(Vector v) const {
        Vector result{};
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)result[i]+=forward[j][i]*v[j];
        return result;
    }
};
EquatorialFrame equatorial_frame(Instant<TT> t,const earth::Nutation& n);

}  // namespace spp::detail
