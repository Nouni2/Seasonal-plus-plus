// Full ELP/MPP02 evaluator. Parameters, arguments and precession equations
// adapted from Yuk Tung Liu's GPL-3.0 implementation (see data/elp_mpp02).
// Series storage and summation use immutable compiled data and no runtime I/O.
#include "spp/model/elp_mpp02.hpp"

#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <type_traits>

#include "../core/split.hpp"
#include "elp_mpp02/arguments.hpp"
#include "elp_mpp02/terms.hpp"
#include "spp/elp_mpp02_terms.inc"

namespace spp::model::elp_mpp02 {
namespace {

using namespace detail;

struct ValueRate { double value = 0.0; double rate = 0.0; };

template<bool WithVelocity>
ValueRate main_sum(std::span<const MainTerm> terms, const std::array<double, 13>& args,
                   [[maybe_unused]] const std::array<double, 13>& rates, const Factors& factors, bool distance) {
    const std::array<double, 6> adjustment = {
        distance ? factors.fA : 1.0, factors.fB1, factors.fB2,
        factors.fB3, factors.fB4, factors.fB5
    };
    ValueRate result;
    for (const MainTerm& term : terms) {
        double amplitude = term.coefficient[0] * adjustment[0];
        for (std::size_t i = 1; i < adjustment.size(); ++i)
            amplitude += term.coefficient[i] * adjustment[i];
        double phase = 0.0;
        [[maybe_unused]] double rate = 0.0;
        for (std::size_t i = 0; i < term.argument.size(); ++i) {
            phase += term.argument[i] * args[i];
            if constexpr (WithVelocity) rate += term.argument[i] * rates[i];
        }
        if (distance) {
            result.value += amplitude * std::cos(phase);
            if constexpr (WithVelocity) result.rate -= amplitude * std::sin(phase) * rate;
        } else {
            result.value += amplitude * std::sin(phase);
            if constexpr (WithVelocity) result.rate += amplitude * std::cos(phase) * rate;
        }
    }
    return result;
}

template<bool WithVelocity>
ValueRate perturbation_sum(std::span<const PerturbationTerm> terms,
                          const std::array<double, 13>& args,
                          [[maybe_unused]] const std::array<double, 13>& rates) {
    ValueRate result;
    for (const PerturbationTerm& term : terms) {
        double phase = term.phase_rad;
        [[maybe_unused]] double rate = 0.0;
        for (std::size_t i = 0; i < term.argument.size(); ++i) {
            phase += term.argument[i] * args[i];
            if constexpr (WithVelocity) rate += term.argument[i] * rates[i];
        }
        result.value += term.amplitude * std::sin(phase);
        if constexpr (WithVelocity) result.rate += term.amplitude * std::cos(phase) * rate;
    }
    return result;
}

std::array<double, 13> components(const Arguments& args) {
    return {args.D, args.F, args.L, args.Lp, args.Me, args.Ve, args.EM,
            args.Ma, args.Ju, args.Sa, args.Ur, args.Ne, args.zeta};
}

template<bool WithVelocity, std::size_t N>
ValueRate perturbed(ValueRate main, const std::array<std::span<const PerturbationTerm>, N>& series,
                    const std::array<double, 13>& args, const std::array<double, 13>& rates,
                    double t, bool angular) {
    double power = 1.0;
    [[maybe_unused]] double previous_power = 1.0;
    for (std::size_t i = 0; i < N; ++i) {
        const ValueRate term = perturbation_sum<WithVelocity>(series[i], args, rates);
        const double value = term.value * power;
        main.value += angular && i > 0 ? mod2pi(value) : value;
        if constexpr (WithVelocity) {
            main.rate += term.rate * power;
            if (i > 0)
                main.rate += static_cast<double>(i) * term.value * previous_power / days_per_century;
        }
        previous_power = power;
        power *= t;
    }
    return main;
}

template<bool WithVelocity>
using Result = std::conditional_t<WithVelocity, GeocentricEclipticJ2000State, GeocentricEclipticJ2000Position>;

template<bool WithVelocity>
Result<WithVelocity> evaluate(double t, Fit fit) {
    Parameters parameters{};
    Factors factors{};
    setup_parameters(fit == Fit::LLR ? 0 : 1, parameters, factors);
    Arguments arguments{}, derivatives{};
    compute_Elp_arguments(t, parameters, arguments);
    if constexpr (WithVelocity) compute_Elp_arguments_dot(t, parameters, derivatives);
    const auto args = components(arguments);
    const auto rates = components(derivatives);
    const std::array<std::span<const PerturbationTerm>, 4> longitude = {
        data::pert_long0, data::pert_long1, data::pert_long2, data::pert_long3
    };
    const std::array<std::span<const PerturbationTerm>, 3> latitude = {
        data::pert_lat0, data::pert_lat1, data::pert_lat2
    };
    const std::array<std::span<const PerturbationTerm>, 4> distance = {
        data::pert_dist0, data::pert_dist1, data::pert_dist2, data::pert_dist3
    };
    auto lon = main_sum<WithVelocity>(data::main_long, args, rates, factors, false);
    lon.value += arguments.W1;
    if constexpr (WithVelocity) lon.rate += derivatives.W1;
    lon = perturbed<WithVelocity>(lon, longitude, args, rates, t, true);
    const auto lat = perturbed<WithVelocity>(main_sum<WithVelocity>(data::main_lat, args, rates, factors, false),
                               latitude, args, rates, t, true);
    auto radius = perturbed<WithVelocity>(main_sum<WithVelocity>(data::main_dist, args, rates, factors, true),
                             distance, args, rates, t, false);
    constexpr double radius_scale = 384747.961370173 / 384747.980674318;
    radius.value *= radius_scale;
    if constexpr (WithVelocity) radius.rate *= radius_scale;

    const double cl = std::cos(lon.value), sl = std::sin(lon.value);
    const double cb = std::cos(lat.value), sb = std::sin(lat.value);
    const double r = radius.value;
    const std::array<double, 3> p = {r*cl*cb, r*sl*cb, r*sb};

    // ELP's rotation from moving ecliptic to mean dynamical J2000 ecliptic.
    // Include its derivative in the velocity, not just the rotated rate.
    const double t2 = t*t, t3 = t*t2, t4 = t2*t2, t5 = t2*t3;
    const double P = 0.10180391e-4*t + 0.47020439e-6*t2 - 0.5417367e-9*t3
                   - 0.2507948e-11*t4 + 0.463486e-14*t5;
    const double Q = -0.113469002e-3*t + 0.12372674e-6*t2 + 0.12654170e-8*t3
                   - 0.1371808e-11*t4 - 0.320334e-14*t5;
    const double s = std::sqrt(1.0 - P*P - Q*Q);
    const std::array<std::array<double, 3>, 3> rotation = {{
        {1.0 - 2*P*P, 2*P*Q, 2*P*s},
        {2*P*Q, 1.0 - 2*Q*Q, -2*Q*s},
        {-2*P*s, 2*Q*s, 1.0 - 2*P*P - 2*Q*Q}
    }};
    Result<WithVelocity> result{};
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            result.position_km[i] += rotation[i][j] * p[j];
    if constexpr (WithVelocity) {
        const double rd = radius.rate;
        const std::array<double, 3> v = {
            rd*cl*cb - r*sl*cb*lon.rate - r*cl*sb*lat.rate,
            rd*sl*cb + r*cl*cb*lon.rate - r*sl*sb*lat.rate,
            rd*sb + r*cb*lat.rate
        };
        const double Pd = (0.10180391e-4 + 0.94040878e-6*t - 1.6252101e-9*t2
                        - 1.0031792e-11*t3 + 2.31743e-14*t4) / days_per_century;
        const double Qd = (-0.113469002e-3 + 0.24745348e-6*t + 0.3796251e-8*t2
                        - 0.5487232e-11*t3 - 1.60167e-14*t4) / days_per_century;
        const double sd = -(P*Pd + Q*Qd) / s;
        const std::array<std::array<double, 3>, 3> rotation_rate = {{
            {-4*P*Pd, 2*(Pd*Q + P*Qd), 2*(Pd*s + P*sd)},
            {2*(Pd*Q + P*Qd), -4*Q*Qd, -2*(Qd*s + Q*sd)},
            {-2*(Pd*s + P*sd), 2*(Qd*s + Q*sd), -4*P*Pd - 4*Q*Qd}
        }};
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                result.velocity_km_per_day[i] += rotation[i][j] * v[j] + rotation_rate[i][j] * p[j];
    }
    return result;
}

Instant<TDB> checked_epoch(Instant<TDB> t) {
    if (!std::isfinite(t.frac))
        throw std::invalid_argument("ELP/MPP02: TDB fraction must be finite");
    t = spp::detail::normalized(t);
    // Check the split before subtracting the epoch, including extreme int64
    // days and fractions that cannot safely be normalized.
    constexpr std::int64_t first = Ephemeris::first_epoch.day;
    constexpr std::int64_t last = Ephemeris::last_epoch.day;
    if (t.frac < 0.0 || t.frac >= 1.0 || t.day < first || t.day > last
        || (t.day == last && t.frac > 0.0))
        throw std::out_of_range("ELP/MPP02: epoch outside -31..11 Julian centuries from J2000");
    return t;
}

}  // namespace

Ephemeris::Ephemeris(Fit fit) : fit_(fit) {
    if (fit != Fit::DE405 && fit != Fit::LLR)
        throw std::invalid_argument("ELP/MPP02: unknown fitted parameter set");
}

GeocentricEclipticJ2000State Ephemeris::operator()(Instant<TDB> t) const {
    return evaluate<true>(checked_epoch(t).centuries_since_j2000(), fit_);
}

GeocentricEclipticJ2000Position Ephemeris::position(Instant<TDB> t) const {
    return evaluate<false>(checked_epoch(t).centuries_since_j2000(), fit_);
}

}  // namespace spp::model::elp_mpp02
