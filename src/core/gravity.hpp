#pragma once

#include "vector.hpp"
#include <array>
#include <stdexcept>

namespace spp::detail {

template<std::size_t Degree>
using GravityCoefficients=std::array<std::array<double,Degree+1>,Degree+1>;

// Closed degree-2 field, also used for distant point-mass/figure reactions.
inline Vector quadrupole_acceleration(Vector r,double radius,double mu,
    double c20,double c21,double s21,double c22,double s22) {
    const Vector q{(-0.5*c20+3*c22)*r[0]+3*s22*r[1]+1.5*c21*r[2],
        3*s22*r[0]+(-0.5*c20-3*c22)*r[1]+1.5*s21*r[2],
        1.5*c21*r[0]+1.5*s21*r[1]+c20*r[2]};
    const double r2=dot(r,r);
    return (mu*radius*radius/(r2*r2*std::sqrt(r2)))*(2*q-(5*dot(r,q)/r2)*r);
}

// Fast axisymmetric counterpart, using ordinary Legendre polynomials and
// their derivatives. J_n has the geodetic sign, C_n0=-J_n. pole is a unit vector.
template<std::size_t Size>
Vector zonal_acceleration(Vector r,Vector pole,double radius,double mu,unsigned degree,const std::array<double,Size>& j) {
    if(degree>=Size || degree<2)throw std::invalid_argument("spp gravity: invalid zonal degree");
    const double distance=norm(r),s=dot(r,pole)/distance;
    const Vector radial=(1.0/distance)*r,tangent=pole-s*radial;
    double p0=1,p1=s,d0=0,d1=1,power=radius/distance;
    Vector g{};
    for(unsigned n=2;n<=degree;++n) {
        power*=radius/distance;
        const double p=((2*n-1)*s*p1-(n-1)*p0)/n;
        const double d=((2*n-1)*(p1+s*d1)-(n-1)*d0)/n;
        g=g+(j[n]*power)*((n+1)*p*radial-d*tangent);
        p0=p1;p1=p;d0=d1;d1=d;
    }
    return (mu/(distance*distance))*g;
}

namespace gravity_detail {
using spp::detail::operator+;
using spp::detail::operator-;
using spp::detail::operator*;
struct Term {double value=0;Vector gradient{};};
inline Term operator+(const Term& a,const Term& b){return {a.value+b.value,a.gradient+b.gradient};}
inline Term operator-(const Term& a,const Term& b){return {a.value-b.value,a.gradient-b.gradient};}
inline Term operator*(double s,const Term& a){return {s*a.value,s*a.gradient};}
inline Term operator*(const Term& a,const Term& b){return {a.value*b.value,a.value*b.gradient+b.value*a.gradient};}
}

// Gradient of mu/R sum (C_nm V_nm + S_nm W_nm), excluding degrees 0/1.
// Unnormalized real solid harmonics, no Condon-Shortley sign. Cartesian
// recurrences and their analytic gradients are regular at both poles.
// r and R use the same distance unit; mu fixes acceleration units.
template<std::size_t MaximumDegree>
Vector harmonic_acceleration(Vector r,double radius,double mu,unsigned degree,
    const GravityCoefficients<MaximumDegree>& c,const GravityCoefficients<MaximumDegree>& s) {
    using namespace gravity_detail;
    if(degree>MaximumDegree || degree<2 || !std::isfinite(radius) || radius<=0 || !std::isfinite(mu))
        throw std::invalid_argument("spp gravity: invalid degree, radius or mass parameter");
    const double distance=norm(r),r2=dot(r,r);
    if(!std::isfinite(distance) || distance<=0 || !std::isfinite(r2))throw std::invalid_argument("spp gravity: invalid position");
    const double inverse=1/r2,k=radius*inverse;
    std::array<Term,3> coordinate{};
    for(std::size_t j=0;j<3;++j) {
        coordinate[j].value=k*r[j];coordinate[j].gradient=(-2*k*r[j]*inverse)*r;
        coordinate[j].gradient[j]+=k;
    }
    const auto& x=coordinate[0];const auto& y=coordinate[1];const auto& z=coordinate[2];
    const Term rho2{radius*radius*inverse,(-2*radius*radius*inverse*inverse)*r};
    std::array<std::array<Term,MaximumDegree+1>,MaximumDegree+1> v{},w{};
    v[0][0]={radius/distance,(-radius/(distance*r2))*r};
    Term potential;
    for(unsigned n=1;n<=degree;++n) {
        const double diagonal=2.0*n-1;
        v[n][n]=diagonal*(x*v[n-1][n-1]-y*w[n-1][n-1]);
        w[n][n]=diagonal*(x*w[n-1][n-1]+y*v[n-1][n-1]);
        for(unsigned m=0;m<n;++m) {
            const double a=(2.0*n-1)/(n-m),b=(n+m-1.0)/(n-m);
            v[n][m]=a*(z*v[n-1][m]);w[n][m]=a*(z*w[n-1][m]);
            if(m+2<=n) {
                v[n][m]=v[n][m]-b*(rho2*v[n-2][m]);
                w[n][m]=w[n][m]-b*(rho2*w[n-2][m]);
            }
        }
        if(n>=2)for(unsigned m=0;m<=n;++m)potential=potential+c[n][m]*v[n][m]+s[n][m]*w[n][m];
    }
    return (mu/radius)*potential.gradient;
}

} // namespace spp::detail
