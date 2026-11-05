#include "frames.hpp"

#include "spp/core/constants.hpp"
#include "spp/core/rotation.hpp"

namespace spp::detail {
namespace {

// Fukushima-Williams P03 angles (Hilton et al. 2006). Coefficients checked
// against ERFA pfw06/obl06; see data/licenses/ERFA.txt for attribution.
struct Angles { double gamma, phi, psi, eps; };
Angles angles(double t) {
    constexpr double as2r = deg2rad / 3600.0;
    return {
        (-0.052928 + t*(10.556378 + t*(0.4932044 + t*(-0.00031238 + t*(-0.000002788 + t*0.0000000260))))) * as2r,
        (84381.412819 + t*(-46.811016 + t*(0.0511268 + t*(0.00053289 + t*(-0.000000440 - t*0.0000000176))))) * as2r,
        (-0.041775 + t*(5038.481484 + t*(1.5584175 + t*(-0.00018522 + t*(-0.000026452 - t*0.0000000148))))) * as2r,
        (84381.406 + t*(-46.836769 + t*(-0.0001831 + t*(0.00200340 + t*(-0.000000576 - t*0.0000000434))))) * as2r
    };
}

Vector precess(Vector v, Angles a) {
    return rotate_x(rotate_z(rotate_x(rotate_z(v, -a.gamma), -a.phi), a.psi), a.eps);
}
Vector unprecess(Vector v, Angles a) {
    return rotate_z(rotate_x(rotate_z(rotate_x(v, -a.eps), -a.psi), a.phi), a.gamma);
}
constexpr double eps_j2000 = 84381.448 * deg2rad / 3600.0;

}  // namespace

Vector ecliptic_j2000_to_mean_equatorial(Vector v, Instant<TT> t) {
    // Remove the J2000 bias: P(t) = B(t) transpose(B(0)).
    return precess(unprecess(rotate_x(v, eps_j2000), angles(0.0)), angles(t.centuries_since_j2000()));
}
Vector mean_equatorial_to_ecliptic_j2000(Vector v, Instant<TT> t) {
    return rotate_x(precess(unprecess(v, angles(t.centuries_since_j2000())), angles(0.0)), -eps_j2000);
}
Vector ecliptic_j2000_to_equatorial(Vector v, Instant<TT> t, const earth::Nutation& n) {
    v = ecliptic_j2000_to_mean_equatorial(v, t);
    return rotate_x(rotate_z(rotate_x(v, -n.mean_obliquity_deg * deg2rad), n.dpsi_deg * deg2rad),
                    n.true_obliquity_deg * deg2rad);
}
Vector equatorial_to_ecliptic_j2000(Vector v, Instant<TT> t, const earth::Nutation& n) {
    v = rotate_x(rotate_z(rotate_x(v, -n.true_obliquity_deg * deg2rad), -n.dpsi_deg * deg2rad),
                 n.mean_obliquity_deg * deg2rad);
    return mean_equatorial_to_ecliptic_j2000(v, t);
}

EquatorialFrame equatorial_frame(Instant<TT> t,const earth::Nutation& n) {
    EquatorialFrame result{{{{1,0,0},{0,1,0},{0,0,1}}}};
    auto rx=[&](double angle) {
        const double c=std::cos(angle),s=std::sin(angle);
        for(std::size_t j=0;j<3;++j) {
            const double y=result.forward[1][j],z=result.forward[2][j];
            result.forward[1][j]=c*y-s*z;result.forward[2][j]=s*y+c*z;
        }
    };
    auto rz=[&](double angle) {
        const double c=std::cos(angle),s=std::sin(angle);
        for(std::size_t j=0;j<3;++j) {
            const double x=result.forward[0][j],y=result.forward[1][j];
            result.forward[0][j]=c*x-s*y;result.forward[1][j]=s*x+c*y;
        }
    };
    const auto zero=angles(0),date=angles(t.centuries_since_j2000());
    rx(eps_j2000);
    rx(-zero.eps);rz(-zero.psi);rx(zero.phi);rz(zero.gamma);
    rz(-date.gamma);rx(-date.phi);rz(date.psi);rx(date.eps);
    rx(-n.mean_obliquity_deg*deg2rad);rz(n.dpsi_deg*deg2rad);rx(n.true_obliquity_deg*deg2rad);
    return result;
}

}  // namespace spp::detail

namespace spp::frames {
FrameRotation<frame::MeanEclipticJ2000,frame::Icrf> mean_ecliptic_j2000_to_icrf() {
    // At J2000, unprecess is the inverse of the ICRF -> dynamical J2000
    // frame-bias rotation. The existing date transformations explicitly
    // cancel this bias, so their mean-of-date operation stays unchanged.
    FrameRotation<frame::MeanEclipticJ2000,frame::Icrf> result;
    for(std::size_t j=0;j<3;++j) {
        detail::Vector axis{};axis[j]=1;
        const auto v=detail::unprecess(detail::rotate_x(axis,detail::eps_j2000),detail::angles(0));
        for(std::size_t i=0;i<3;++i)result.matrix[i][j]=v[i];
    }
    return result;
}
}  // namespace spp::frames
