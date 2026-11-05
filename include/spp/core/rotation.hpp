#pragma once

#include "spp/core/state.hpp"

namespace spp {

using RotationMatrix=std::array<std::array<double,3>,3>;

// Rotation of axes, preserving the origin tag. All quantities must refer to
// the same epoch. The derivative is per TDB day, as are state velocities.
template<class From,class To>
struct FrameRotation {
    RotationMatrix matrix{};
    RotationMatrix derivative_per_day{};

    template<class Origin>
    CartesianPosition<Origin,To> operator()(const CartesianPosition<Origin,From>& p) const {
        CartesianPosition<Origin,To> result{};
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
            result.position_km[i]+=matrix[i][j]*p.position_km[j];
        return result;
    }
    template<class Origin>
    CartesianState<Origin,To> operator()(const CartesianState<Origin,From>& s) const {
        CartesianState<Origin,To> result{};
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
            result.position_km[i]+=matrix[i][j]*s.position_km[j];
            result.velocity_km_per_day[i]+=matrix[i][j]*s.velocity_km_per_day[j]
                +derivative_per_day[i][j]*s.position_km[j];
        }
        return result;
    }
    FrameRotation<To,From> inverse() const {
        FrameRotation<To,From> result;
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
            result.matrix[i][j]=matrix[j][i];
            result.derivative_per_day[i][j]=derivative_per_day[j][i];
        }
        return result;
    }
};

template<class From,class Middle,class To>
FrameRotation<From,To> operator*(const FrameRotation<Middle,To>& a,const FrameRotation<From,Middle>& b) {
    FrameRotation<From,To> result;
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)for(std::size_t k=0;k<3;++k) {
        result.matrix[i][j]+=a.matrix[i][k]*b.matrix[k][j];
        result.derivative_per_day[i][j]+=a.derivative_per_day[i][k]*b.matrix[k][j]
            +a.matrix[i][k]*b.derivative_per_day[k][j];
    }
    return result;
}

namespace frames {
// Fixed J2000 tilt (84381.448 arcsec), then removal of the P03 J2000 frame
// bias. The inverse is available through FrameRotation::inverse().
FrameRotation<frame::MeanEclipticJ2000,frame::Icrf> mean_ecliptic_j2000_to_icrf();
}
}  // namespace spp
