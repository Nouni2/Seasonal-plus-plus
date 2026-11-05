#pragma once

#include <concepts>
#include <optional>
#include <vector>

#include "spp/core/state.hpp"
#include "spp/core/rotation.hpp"
#include "spp/core/time.hpp"
#include "spp/core/types.hpp"

namespace spp::moon {

using State = GeocentricEclipticJ2000State;
using Position = GeocentricEclipticJ2000Position;

// An ephemeris implementation supplies the Moon's geometric state in this
// exact frame/origin at TDB. The model may own caches or propagation state;
// the body API imposes no analytical-model-specific initialization scheme.
template<class Model>
concept Ephemeris = requires(Model& model, Instant<TDB> t) {
    { model(t) } -> std::same_as<State>;
};

// Observation does not require a velocity. Models can supply the same typed
// position directly, or retain the existing full-state callable contract.
template<class Model>
concept PositionEphemeris = requires(Model& model, Instant<TDB> t) {
    { model.position(t) } -> std::same_as<Position>;
};
template<class Model>
concept PositionSource = PositionEphemeris<Model> || Ephemeris<Model>;

using BodyRotation=FrameRotation<frame::Icrf,frame::MoonMeanEarth>;
template<class Model>
concept RotationSource = requires(Model& model,Instant<TDB> t) {
    {model(t)} -> std::convertible_to<BodyRotation>;
} || requires(Model& model,Instant<TDB> t) {
    {model(t).rotation} -> std::convertible_to<BodyRotation>;
};

// Default: full ELP/MPP02, DE405/DE406-adjusted parameters, direct evaluation.
State geometric(Instant<TDB> t);
Position position(Instant<TDB> t);

// Explicit model selection, with compile-time dispatch in repeated queries.
// Models with different frames/origins need an explicit adapter.
template<class Model> requires Ephemeris<Model>
State geometric(Instant<TDB> t, Model& model) {
    return model(t);
}

template<class Model> requires PositionSource<Model>
Position position(Instant<TDB> t, Model& model) {
    if constexpr (PositionEphemeris<Model>) return model.position(t);
    else return {model(t).position_km};
}

// Geocentric phase. The longitude cycle (0=new, 90=first quarter, 180=full,
// 270=last quarter) differs from the physical Sun-Moon-Earth phase angle.
struct Phase {
    double cycle_deg;               // apparent ecliptic longitude Moon - Sun, [0,360)
    double elongation_deg;          // angular separation from Sun, [0,180]
    double phase_angle_deg;         // physical illumination angle at reflection: 0=fully illuminated
    double illuminated_fraction;    // [0,1], spherical Moon, excludes eclipses
    bool waxing;                   // cycle < 180 degrees
    const char* name() const;       // crescent/gibbous; "Near ..." within 0.5 deg of a principal phase
};

// Spherical disk as seen from the observer, before differential atmospheric
// distortion. Illumination uses the astrometric outgoing ray and incoming
// apparent Sun at reflection, separately from the geocentric phase cycle.
// Angles are measured eastward from celestial north; subtract
// parallactic_angle_deg to measure from the local zenith. Bright-limb angle
// is undefined at exact alignment; the zenith direction is undefined at the
// zenith/nadir. Flags distinguish those cases from a valid zero-degree angle.
struct Disk {
    double phase_angle_deg;
    double illuminated_fraction;
    double bright_limb_deg;
    double parallactic_angle_deg;
    bool bright_limb_defined;
    bool zenith_defined;
};

// Spherical Mean Earth/Polar Axis surface coordinates at reflection.
// Longitudes are east-positive [-180,180); latitude is planetocentric.
struct Surface {
    double subobserver_lon_deg,subobserver_lat_deg;
    double subsolar_lon_deg,subsolar_lat_deg;
    double pole_position_angle_deg; // east from celestial north
    bool pole_defined;
    // Orthographic display axes: screen right, celestial north, towards viewer.
    // Directions are body-fixed unit vectors; atmospheric distortion omitted.
    RotationMatrix body_to_view;
};

struct Observation {
    Equatorial equatorial;          // apparent topocentric, true equator/equinox of date
    Horizontal horizontal;          // center, refraction according to Atmosphere
    double geometric_alt_deg;       // center, without atmospheric refraction
    double distance_km;             // observer to Moon, with light travel time
    double semidiameter_deg;        // spherical mean radius 1737.4 km
    Phase phase;                    // geocentric, independent of observer/refraction
    Disk disk;                      // topocentric illumination and disk orientation
    Surface surface;                // default IAU2009 orientation, independently replaceable
};

// Observer coordinates without geocentric phase/disk-orientation work.
struct Coordinates {
    Equatorial equatorial;
    Horizontal horizontal;
    double geometric_alt_deg;
    double distance_km;
    double semidiameter_deg;
    double hour_angle_deg;           // [-180,180), positive west of the meridian
};

enum class HorizonKind { Rise, Set, UpperTransit, LowerTransit, Graze };
const char* name(HorizonKind kind);
struct HorizonEvent {
    HorizonKind kind;
    Time time;
    double center_altitude_deg;      // apparent if the configured refraction is enabled
    double azimuth_deg;
};
struct HorizonOptions {
    double altitude_deg=0;           // horizontal obstruction angle; no automatic terrain/dip
    bool upper_limb=true;
    double time_tolerance_s=0.01;    // numerical tolerance, not prediction accuracy
    double contact_tolerance_arcsec=0.001;
};
struct HorizonResult {
    std::vector<HorizonEvent> events; // all events, sorted, in [start,end)
    bool above_at_start=false;       // immediately after start; false for an empty interval
    double duration_hours=0;         // portion of the requested UTC interval above the target
};

enum class PrincipalPhase { NewMoon, FirstQuarter, FullMoon, LastQuarter };
const char* name(PrincipalPhase phase);
struct PhaseEvent {
    PrincipalPhase phase;
    Time time;                      // UTC; same time-scale policy as spp::Time
};

// UTC-bearing Time supplies TT/TDB for the ephemeris and UT1 for rotation.
// Altitude is the center's angle above the horizon, not a physical height.
// Default observation/phase calls share a per-thread, bounded twenty-minute
// ELP/MPP02 Hermite cache. Explicit providers retain their own policies;
// supplying elp_mpp02::Ephemeris requests direct full-series evaluation.
Observation observe(const Time& t, const Observer& obs, const Atmosphere& atm = {});
Surface surface(const Time& t);     // geocentric, default cached orbit and IAU rotation
Surface surface(const Time& t,const Observer& obs);
Coordinates coordinates(const Time& t,const Observer& obs,const Atmosphere& atm = {});
Equatorial apparent(const Time& t);
Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm = {});
Phase phase(const Time& t);
// Principal phases in [start,end), in time order. 0.01 s is a numerical root
// tolerance, not a claim of physical accuracy. Reversed intervals are errors.
std::vector<PhaseEvent> phases(const Time& start, const Time& end, double tolerance_s = 0.01);
// Terrestrial Moon search, including both meridian crossings and stationary
// contacts. Uses a conservative smooth-geometry curvature bound, not the solar
// transit/constant-declination approximation. Refraction is inverted once.
HorizonResult horizon(const Time& start,const Time& end,const Observer& obs,
                      const Atmosphere& atm = {},const HorizonOptions& options = {});

namespace detail {
// A synchronous, non-owning model reference keeps observation algorithms
// independent of ELP, without allocating or copying a provider's caches.
struct ModelRef {
    const void* context;
    Position (*evaluate)(const void*, Instant<TDB>);

    template<class Model> requires PositionSource<Model>
    explicit ModelRef(Model& model) : context(&model), evaluate([](const void* p, Instant<TDB> t) {
        return moon::position(t, *static_cast<Model*>(const_cast<void*>(p)));
    }) {}
    Position operator()(Instant<TDB> t) const { return evaluate(context, t); }
};
struct RotationRef {
    const void* context;
    BodyRotation (*evaluate)(const void*,Instant<TDB>);
    template<class Model> requires RotationSource<Model>
    explicit RotationRef(Model& model):context(&model),evaluate([](const void* p,Instant<TDB> t) {
        auto& source=*static_cast<Model*>(const_cast<void*>(p));
        if constexpr(requires{{source(t)} -> std::convertible_to<BodyRotation>;})return source(t);
        else return source(t).rotation;
    }) {}
    BodyRotation operator()(Instant<TDB> t) const{return evaluate(context,t);}
};
Observation observe(const Time&, const Observer&, const Atmosphere&, ModelRef);
Observation observe(const Time&,const Observer&,const Atmosphere&,ModelRef,RotationRef);
Surface surface(const Time&,const std::optional<Observer>&,ModelRef);
Surface surface(const Time&,const std::optional<Observer>&,ModelRef,RotationRef);
Coordinates coordinates(const Time&,const Observer&,const Atmosphere&,ModelRef);
Horizontal horizontal(const Time&, const Observer&, const Atmosphere&, ModelRef);
Equatorial apparent(const Time&, ModelRef);
Phase phase(const Time&, ModelRef);
std::vector<PhaseEvent> phases(const Time&, const Time&, double tolerance_s, ModelRef);
HorizonResult horizon(const Time&,const Time&,const Observer&,const Atmosphere&,const HorizonOptions&,ModelRef);
}  // namespace detail

template<class Model> requires PositionSource<Model>
Observation observe(const Time& t, const Observer& obs, Model& model, const Atmosphere& atm = {}) {
    return detail::observe(t, obs, atm, detail::ModelRef(model));
}
template<class Model,class Spin> requires PositionSource<Model> && RotationSource<Spin>
Observation observe(const Time& t,const Observer& obs,Model& model,Spin& rotation,const Atmosphere& atm={}) {
    return detail::observe(t,obs,atm,detail::ModelRef(model),detail::RotationRef(rotation));
}
template<class Model> requires PositionSource<Model>
Surface surface(const Time& t,Model& model) {
    return detail::surface(t,std::nullopt,detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
Surface surface(const Time& t,const Observer& obs,Model& model) {
    return detail::surface(t,obs,detail::ModelRef(model));
}
template<class Model,class Spin> requires PositionSource<Model> && RotationSource<Spin>
Surface surface(const Time& t,Model& model,Spin& rotation) {
    return detail::surface(t,std::nullopt,detail::ModelRef(model),detail::RotationRef(rotation));
}
template<class Model,class Spin> requires PositionSource<Model> && RotationSource<Spin>
Surface surface(const Time& t,const Observer& obs,Model& model,Spin& rotation) {
    return detail::surface(t,obs,detail::ModelRef(model),detail::RotationRef(rotation));
}
template<class Model> requires PositionSource<Model>
Coordinates coordinates(const Time& t,const Observer& obs,Model& model,const Atmosphere& atm = {}) {
    return detail::coordinates(t,obs,atm,detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
Equatorial apparent(const Time& t, Model& model) {
    return detail::apparent(t, detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
Horizontal horizontal(const Time& t, const Observer& obs, Model& model, const Atmosphere& atm = {}) {
    return detail::horizontal(t, obs, atm, detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
Phase phase(const Time& t, Model& model) {
    return detail::phase(t, detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
std::vector<PhaseEvent> phases(const Time& start, const Time& end, Model& model, double tolerance_s = 0.01) {
    return detail::phases(start, end, tolerance_s, detail::ModelRef(model));
}
template<class Model> requires PositionSource<Model>
HorizonResult horizon(const Time& start,const Time& end,const Observer& obs,Model& model,
                      const Atmosphere& atm = {},const HorizonOptions& options = {}) {
    return detail::horizon(start,end,obs,atm,options,detail::ModelRef(model));
}

}  // namespace spp::moon
