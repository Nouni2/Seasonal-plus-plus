#pragma once

#include "vector.hpp"
#include <stdexcept>

namespace spp::detail {
// Cubic local state prediction. Its short-lag domain is bounded by the caller;
// no external satellite trajectory or delayed/advanced integration state is used.
inline Vector lagged_position(Vector r,Vector v,Vector a,Vector jerk,double lag_days) {
    return r-lag_days*v+(0.5*lag_days*lag_days)*a-(lag_days*lag_days*lag_days/6)*jerk;
}
inline Vector rotate_about(Vector v,Vector pole,double angle) {
    return std::cos(angle)*v+std::sin(angle)*cross(pole,v)+((1-std::cos(angle))*dot(pole,v))*pole;
}
// One order of a degree-2 tidal potential. The shape tensor is the m=0,1,2
// part of q q^T - |q|^2 I/3 about a unit pole. Summing equal-k zero-lag
// modes reproduces the isotropic equilibrium-tide potential exactly.
// mu includes the inducing mass and any required point-mass reaction factor.
inline Vector tide_acceleration(Vector r,Vector q,Vector pole,unsigned order,double radius,double mu,double love) {
    if(order>2)throw std::invalid_argument("spp tides: order must be 0, 1 or 2");
    const double r2=dot(r,r),q2=dot(q,q),z=dot(r,pole),qz=dot(q,pole);
    if(!std::isfinite(r2) || !std::isfinite(q2) || r2<=0 || q2<=0
        || !std::isfinite(love) || !std::isfinite(mu) || !std::isfinite(radius) || radius<=0)
        throw std::invalid_argument("spp tides: invalid field geometry or parameters");
    const Vector rho=r-z*pole,qp=q-qz*pole;
    Vector shape;
    if(order==0) {
        const double a=qz*qz-q2/3;
        shape=(1.5*a*z)*pole-(0.5*a)*r;
    } else if(order==1)shape=qz*(dot(rho,qp)*pole+z*qp);
    else shape=dot(rho,qp)*qp-(0.5*dot(qp,qp))*rho;
    const double radius2=radius*radius;
    const double scale=1.5*love*mu*(radius2*radius2*radius)/(r2*r2*std::sqrt(r2)*q2*q2*std::sqrt(q2));
    return scale*(2*shape-(5*dot(r,shape)/r2)*r);
}
} // namespace spp::detail
