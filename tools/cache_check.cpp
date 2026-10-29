// cache_check: accuracy of the Chebyshev earth-state cache.
//
// Both the cache and the direct double-precision series are compared against
// the same series evaluated in extended precision (x87 long double, 64-bit
// mantissa) at the exact same instant. Direct double evaluation is itself
// noisy: its VSOP phases B + C t reach ~1e6 rad, where a double resolves only
// ~1e-10 rad. The cache passes when its error is no larger than that of the
// direct evaluation it replaces.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <print>
#include <random>
#include <span>

#include "earth_state.hpp"
#include "model/iau1980_terms.hpp"
#include "spp/model/vsop87d.hpp"
#include "spp/vsop87d_earth_terms.inc"

namespace {

using ld = long double;
namespace vsop = spp::model::vsop87d;

// Where long double is no wider than double (MSVC, ARM64), the reference would
// be no more precise than what it checks.
static_assert(std::numeric_limits<ld>::digits > std::numeric_limits<double>::digits,
              "cache_check needs an extended-precision long double");

constexpr ld kPi = 3.141592653589793238462643383279502884L;

ld series(const std::span<const vsop::Term> (&s)[6], ld t) {
    ld value = 0.0L, t_power = 1.0L;
    for (const auto& terms : s) {
        ld inner = 0.0L;
        for (const auto& term : terms) inner += ld(term.a) * cosl(ld(term.b) + ld(term.c) * t);
        value += inner * t_power;
        t_power *= t;
    }
    return value;
}

ld norm360l(ld x) {
    ld r = fmodl(x, 360.0L);
    return r < 0 ? r + 360.0L : r;
}

// The implemented model (double coefficients), evaluated in extended precision.
struct Truth {
    ld l, b, r, dpsi, deps;
};

Truth truth(std::int64_t days, double frac) {
    const ld T = (ld(days) + ld(frac)) / 36525.0L;
    const ld t = T / 10.0L;
    Truth out{series(vsop::data::L, t), series(vsop::data::B, t), series(vsop::data::R, t), 0, 0};

    const ld T2 = T * T, T3 = T2 * T;
    const ld deg = kPi / 180.0L;
    const ld D  = norm360l(ld(297.85036) + ld(445267.111480) * T - ld(0.0019142) * T2 + T3 / ld(189474.0)) * deg;
    const ld M  = norm360l(ld(357.52772) + ld(35999.050340) * T - ld(0.0001603) * T2 - T3 / ld(300000.0)) * deg;
    const ld Mp = norm360l(ld(134.96298) + ld(477198.867398) * T + ld(0.0086972) * T2 + T3 / ld(56250.0)) * deg;
    const ld F  = norm360l(ld(93.27191) + ld(483202.017538) * T - ld(0.0036825) * T2 + T3 / ld(327270.0)) * deg;
    const ld Om = norm360l(ld(125.04452) - ld(1934.136261) * T + ld(0.0020708) * T2 + T3 / ld(450000.0)) * deg;
    for (const auto& r : spp::model::iau1980::kTerms) {
        const ld arg = r.d * D + r.m * M + r.mp * Mp + r.f * F + r.om * Om;
        out.dpsi += (ld(r.psi_a) + ld(r.psi_at) * T) * sinl(arg);
        out.deps += (ld(r.eps_b) + ld(r.eps_bt) * T) * cosl(arg);
    }
    out.dpsi *= 1e-4L;
    out.deps *= 1e-4L;
    return out;
}

struct Stat {
    const char* name;
    double scale;
    const char* unit;
    double direct = 0.0;
    double cached = 0.0;
    long worse = 0;   // instants where the cache is further from truth than direct

    void add(double c, double d, ld truth) {
        const double ec = double(fabsl(ld(c) - truth));
        const double ed = double(fabsl(ld(d) - truth));
        direct = std::max(direct, ed);
        cached = std::max(cached, ec);
        if (ec > ed) ++worse;
    }
};

// The design choices of the cache, measured: the same Chebyshev fit as
// src/earth_state.cpp, rebuilt with a different segment length or without the
// offset, evaluated in double, against the extended-precision series.
struct Variant {
    const char* name;
    std::int64_t segment_days;
    bool offset;
    ld max_err[3] = {0, 0, 0};   // L, B (rad), R (AU)
};

void evaluate_variant(Variant& v, std::int64_t days, double frac, const Truth& x) {
    constexpr int kDegree = 24, kNodes = kDegree + 1;
    const std::int64_t index = (days >= 0 ? days : days - v.segment_days + 1) / v.segment_days;
    const double start = double(index * v.segment_days), half = 0.5 * double(v.segment_days);
    double f[3][kNodes];
    for (int k = 0; k < kNodes; ++k) {
        const double xk = std::cos(double(kPi) * (k + 0.5) / kNodes);
        const auto e = spp::detail::earth_state_direct((start + half * (xk + 1.0)) / 36525.0);
        f[0][k] = e.l_rad;
        f[1][k] = e.b_rad;
        f[2][k] = e.r_au;
    }
    const double local = double(days - index * v.segment_days) + frac;
    const double xx = local * (2.0 / double(v.segment_days)) - 1.0;
    const ld truth_q[3] = {x.l, x.b, x.r};
    for (int q = 0; q < 3; ++q) {
        const double off = v.offset ? f[q][kNodes / 2] : 0.0;
        double c[kNodes];
        for (int j = 0; j < kNodes; ++j) {
            double sum = 0.0;
            for (int k = 0; k < kNodes; ++k) sum += (f[q][k] - off) * std::cos(double(kPi) * j * (k + 0.5) / kNodes);
            c[j] = (2.0 / kNodes) * sum;
        }
        c[0] *= 0.5;
        double b1 = 0.0, b2 = 0.0;
        for (int j = kDegree; j >= 1; --j) {
            const double b0 = 2.0 * xx * b1 - b2 + c[j];
            b2 = b1;
            b1 = b0;
        }
        const double value = off + (xx * b1 - b2 + c[0]);
        v.max_err[q] = std::max(v.max_err[q], fabsl(ld(value) - truth_q[q]));
    }
}

}  // namespace

int main() {
    constexpr double arcsec_per_rad = 206264.80624709636;
    Stat stats[] = {
        {"L", arcsec_per_rad, "\""}, {"B", arcsec_per_rad, "\""}, {"R", 1.0, "AU"},
        {"dpsi", 1.0, "\""},         {"deps", 1.0, "\""},
    };

    std::mt19937_64 rng(2026);
    const std::int64_t first = 1355808 - spp::jd_j2000_day;   // 1000 BC
    const std::int64_t last = 2816788 - spp::jd_j2000_day;    // 3000 AD
    std::uniform_int_distribution<std::int64_t> day_dist(first, last);
    std::uniform_real_distribution<double> frac_dist(0.0, 1.0);

    long n = 0;
    auto check = [&](std::int64_t days, double frac) {
        const spp::Instant<spp::TT> t{spp::jd_j2000_day + days, frac};
        const auto c = spp::detail::earth_state(t);
        const auto d = spp::detail::earth_state_direct(t.centuries_since_j2000());
        const Truth x = truth(days, frac);
        stats[0].add(c.l_rad, d.l_rad, x.l);
        stats[1].add(c.b_rad, d.b_rad, x.b);
        stats[2].add(c.r_au, d.r_au, x.r);
        stats[3].add(c.dpsi_arcsec, d.dpsi_arcsec, x.dpsi);
        stats[4].add(c.deps_arcsec, d.deps_arcsec, x.deps);
        ++n;
    };

    // Random instants over the whole range, dense sampling inside random
    // segments, and both sides of segment boundaries (16-day segments).
    for (int i = 0; i < 10000; ++i) check(day_dist(rng), frac_dist(rng));
    for (int s = 0; s < 100; ++s) {
        const std::int64_t base = day_dist(rng) & ~std::int64_t{15};
        for (int k = 0; k < 100; ++k) check(base + k % 16, frac_dist(rng));
    }
    for (std::int64_t d = -800; d <= 800; d += 16) {
        check(d, 0.0);
        check(d - 1, std::nextafter(1.0, 0.0));
    }

    std::println("{} instants, 1000 BC .. 3000 AD; error vs extended-precision evaluation", n);
    std::println("{:<6}{:>18}{:>18}{:>16}", "", "direct (double)", "cached", "cache worse at");
    // Pass: the cache is no further from the truth than direct evaluation, up
    // to rounding noise (a factor of 2 on the worst case).
    bool pass = true;
    for (const Stat& s : stats) {
        const bool ok = s.cached <= 2.0 * s.direct;
        pass = pass && ok;
        std::println("{:<6}{:>13.3e} {:<4}{:>13.3e} {:<4}{:>10.1f} %{}", s.name, s.direct * s.scale, s.unit,
                     s.cached * s.scale, s.unit, 100.0 * double(s.worse) / double(n), ok ? "" : "   FAIL");
    }

    // Design choices, for the README: the library's fit, the fit without the
    // offset, and 32-day segments, at 2,000 random instants.
    Variant variants[] = {{"16 d, offset (library)", 16, true}, {"16 d, no offset", 16, false},
                          {"32 d, offset", 32, true}};
    std::mt19937_64 rng2(77);
    for (int i = 0; i < 2000; ++i) {
        const std::int64_t days = day_dist(rng2);
        const double frac = frac_dist(rng2);
        const Truth x = truth(days, frac);
        for (Variant& v : variants) evaluate_variant(v, days, frac, x);
    }
    std::println("\nfit variants, 2000 random instants; max error vs extended-precision evaluation");
    std::println("{:<26}{:>14}{:>14}{:>14}", "", "L", "B", "R");
    for (const Variant& v : variants)
        std::println("{:<26}{:>12.3e} \"{:>12.3e} \"{:>11.3e} AU", v.name, double(v.max_err[0]) * arcsec_per_rad,
                     double(v.max_err[1]) * arcsec_per_rad, double(v.max_err[2]));
    std::println("no offset: L error x{:.2f} ({:.1f} bits)", double(variants[1].max_err[0] / variants[0].max_err[0]),
                 std::log2(double(variants[1].max_err[0] / variants[0].max_err[0])));

    std::println("{}", pass ? "PASS" : "FAIL");
    return pass ? 0 : 1;
}
