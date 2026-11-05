#pragma once

#include "spp/core/state.hpp"
#include "spp/core/time.hpp"
#include "spp/model/hermite.hpp"

namespace spp::model::elp_mpp02 {

// Chapront & Francou (2003), ELP/MPP02. Both published parameter adjustments
// are supported; neither is a fit to DE440/DE441.
enum class Fit { DE405, LLR };

// Full, untruncated analytical series. Immutable configuration; no runtime
// data files, propagation, interpolation cache, allocation or shared state.
class Ephemeris {
public:
    explicit Ephemeris(Fit fit = Fit::DE405);

    static constexpr Instant<TDB> first_epoch{1319270, 0.0};
    static constexpr Instant<TDB> last_epoch{2853320, 0.0};

    // Geometric Moon relative to Earth's center, in the mean dynamical
    // ecliptic and equinox of J2000. Position: km; velocity: km/TDB day.
    // t must be TDB, not TT relabeled as TDB. No apparent-place corrections.
    // Throws std::invalid_argument for invalid fit/non-finite fractions and
    // std::out_of_range outside -31 <= T <= 11 Julian centuries from J2000.
    // This computational interval is not a uniform astronomical accuracy claim.
    GeocentricEclipticJ2000State operator()(Instant<TDB> t) const;
    // Same position, coefficients, frame and input/error contract as the full
    // state. Skips all series derivatives and the velocity transformation.
    GeocentricEclipticJ2000Position position(Instant<TDB> t) const;
    Fit fit() const { return fit_; }

private:
    Fit fit_;
};

// Twenty-minute cubic segments sampled from the full analytical state and its
// velocity. Configuration and the bounded cache belong to this instance;
// different fits/providers cannot share stale entries. Not thread safe: keep
// one instance per worker. Direct Ephemeris remains available for references.
class CachedEphemeris : public HermiteEphemeris<Ephemeris,TDB,72> {
public:
    explicit CachedEphemeris(Fit fit = Fit::DE405)
        : HermiteEphemeris(Ephemeris(fit), Ephemeris::first_epoch, Ephemeris::last_epoch) {}
    Fit fit() const { return source().fit(); }
};

}  // namespace spp::model::elp_mpp02
