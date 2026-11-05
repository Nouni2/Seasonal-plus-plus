#pragma once

#include <array>
#include <cmath>

namespace spp::detail {

using Vector = std::array<double, 3>;

inline Vector operator+(const Vector& a, const Vector& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
inline Vector operator-(const Vector& a, const Vector& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Vector operator*(double s, const Vector& v) { return {s * v[0], s * v[1], s * v[2]}; }
inline double dot(const Vector& a, const Vector& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline Vector cross(const Vector& a, const Vector& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
inline double norm(const Vector& v) { return std::hypot(v[0], v[1], v[2]); }
inline double angle(const Vector& a, const Vector& b) { return std::atan2(norm(cross(a,b)),dot(a,b)); }
inline Vector rotate_x(Vector v, double angle) {
    const double c = std::cos(angle), s = std::sin(angle);
    return {v[0], c*v[1] - s*v[2], s*v[1] + c*v[2]};
}
inline Vector rotate_z(Vector v, double angle) {
    const double c = std::cos(angle), s = std::sin(angle);
    return {c*v[0] - s*v[1], s*v[0] + c*v[1], v[2]};
}

}  // namespace spp::detail
