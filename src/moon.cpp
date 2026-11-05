#include "spp/moon.hpp"

#include <cmath>

#include "spp/core/angles.hpp"
#include "spp/model/elp_mpp02.hpp"

namespace spp::moon {
namespace {
model::elp_mpp02::CachedEphemeris& observation_model() {
    thread_local model::elp_mpp02::CachedEphemeris model;
    return model;
}
}  // namespace

State geometric(Instant<TDB> t) {
    const model::elp_mpp02::Ephemeris ephemeris;
    return geometric(t, ephemeris);
}

Position position(Instant<TDB> t) {
    const model::elp_mpp02::Ephemeris model;
    return position(t, model);
}

Observation observe(const Time& t, const Observer& obs, const Atmosphere& atm) {
    return observe(t, obs, observation_model(), atm);
}
Surface surface(const Time& t) {return surface(t,observation_model());}
Surface surface(const Time& t,const Observer& obs) {return surface(t,obs,observation_model());}
Coordinates coordinates(const Time& t,const Observer& obs,const Atmosphere& atm) {
    return coordinates(t,obs,observation_model(),atm);
}
HorizonResult horizon(const Time& start,const Time& end,const Observer& obs,
                      const Atmosphere& atm,const HorizonOptions& options) {
    return horizon(start,end,obs,observation_model(),atm,options);
}
Equatorial apparent(const Time& t) {
    return apparent(t, observation_model());
}
Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm) {
    return horizontal(t, obs, observation_model(), atm);
}
Phase phase(const Time& t) {
    return phase(t, observation_model());
}
std::vector<PhaseEvent> phases(const Time& start, const Time& end, double tolerance_s) {
    return phases(start, end, observation_model(), tolerance_s);
}

const char* Phase::name() const {
    if (!std::isfinite(cycle_deg)) return "Unknown";
    const double c = norm360(cycle_deg);
    if (c < 0.5 || c > 359.5) return "Near New Moon";
    if (std::abs(c - 90.0) < 0.5) return "Near First Quarter";
    if (std::abs(c - 180.0) < 0.5) return "Near Full Moon";
    if (std::abs(c - 270.0) < 0.5) return "Near Last Quarter";
    if (c < 90.0) return "Waxing Crescent";
    if (c < 180.0) return "Waxing Gibbous";
    if (c < 270.0) return "Waning Gibbous";
    return "Waning Crescent";
}

const char* name(PrincipalPhase p) {
    switch (p) {
        case PrincipalPhase::NewMoon: return "New Moon";
        case PrincipalPhase::FirstQuarter: return "First Quarter";
        case PrincipalPhase::FullMoon: return "Full Moon";
        case PrincipalPhase::LastQuarter: return "Last Quarter";
    }
    return "Unknown";
}

const char* name(HorizonKind kind) {
    switch(kind) {
    case HorizonKind::Rise:return "Moonrise";
    case HorizonKind::Set:return "Moonset";
    case HorizonKind::UpperTransit:return "Upper transit";
    case HorizonKind::LowerTransit:return "Lower transit";
    case HorizonKind::Graze:return "Grazing contact";
    }
    return "Unknown";
}

}  // namespace spp::moon
