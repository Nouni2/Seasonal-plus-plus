#pragma once

#include <cmath>
#include <stdexcept>
#include "spp/core/types.hpp"

namespace spp::moon::detail {
inline constexpr double radius_km=1737.4;

inline void validate(const Observer& obs,const Atmosphere& atm) {
    if(!std::isfinite(obs.lat_deg) || std::abs(obs.lat_deg)>90 || !std::isfinite(obs.lon_deg)
        || !std::isfinite(obs.elev_m) || obs.elev_m< -12000 || obs.elev_m>1e7)
        throw std::invalid_argument("spp::moon: invalid terrestrial observer (latitude, longitude or elevation)");
    if(!std::isfinite(atm.temp_c) || atm.temp_c<=-273
        || (atm.pressure_mbar && (!std::isfinite(*atm.pressure_mbar) || *atm.pressure_mbar<0)))
        throw std::invalid_argument("spp::moon: invalid atmosphere");
}
}  // namespace spp::moon::detail
