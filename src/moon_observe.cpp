#include "spp/moon.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/aberration.hpp"
#include "core/frames.hpp"
#include "core/root_refine.hpp"
#include "core/split.hpp"
#include "moon_internal.hpp"
#include "spp/core/angles.hpp"
#include "spp/earth.hpp"
#include "spp/sky.hpp"
#include "spp/sun.hpp"
#include "spp/model/iau_moon.hpp"

namespace spp::moon::detail {
namespace {

using namespace spp::detail;
constexpr double c_km_day = speed_of_light_m_s * seconds_per_day / 1000.0;
constexpr double moon_radius_km = radius_km;

template<class Scale>
Instant<Scale> shifted(Instant<Scale> t, double days) {
    t.frac += days;
    return normalized(t);
}

Vector spherical(double lon_deg, double lat_deg, double r) {
    const double l = lon_deg * deg2rad, b = lat_deg * deg2rad;
    return {r * std::cos(b) * std::cos(l), r * std::cos(b) * std::sin(l), r * std::sin(b)};
}
Equatorial equatorial(Vector v) {
    return {norm360(std::atan2(v[1], v[0]) * rad2deg),
            std::atan2(v[2], std::hypot(v[0], v[1])) * rad2deg, norm(v) / astronomical_unit_km};
}

Vector earth_position(Instant<TT> tt) {
    // The existing geometric Sun is FK5 mean ecliptic/equinox OF DATE.
    // Reflect it, then undo the date-dependent rotation before differencing.
    const auto sun = sun::geometric(tt);
    const auto v = spherical(sun.lon_deg, sun.lat_deg, -sun.dist_au * astronomical_unit_km);
    return mean_equatorial_to_ecliptic_j2000(rotate_x(v, earth::mean_obliquity_deg(tt) * deg2rad), tt);
}

Vector model_position(ModelRef model, Instant<TDB> t) {
    const auto p = model(t).position_km;
    if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]))
        throw std::invalid_argument("spp::moon: ephemeris returned a non-finite position");
    return p;
}

struct Context {
    Instant<TT> tt;
    Instant<TDB> tdb;
    Vector moon_position_km;
    earth::Nutation nutation;
    EquatorialFrame frame;
    Vector earth_position_km;
    Vector earth_velocity_km_day;

    // Check the provider's epoch/position before evaluating Earth orientation.
    // The reception position seeds both rays in a complete observation.
    Context(const Time& t, ModelRef model) : tt(t.tt()), tdb(tdb_from_tt(tt)),
        moon_position_km(model_position(model, tdb)), nutation(earth::nutation(tt)),frame(equatorial_frame(tt,nutation)),
        earth_position_km(earth_position(tt)) {
        constexpr double step = 1.0 / 128.0;
        earth_velocity_km_day = (0.5 / step) * (earth_position(shifted(tt, step)) - earth_position(shifted(tt, -step)));
    }
};

struct Ray {
    Vector apparent_j2000;
    Vector astrometric_j2000;
    Vector heliocentric_at_emission;
    Vector moon_velocity_km_day;
    Vector observer_velocity_km_day;
    double distance_km;
    Instant<TDB> emission_tdb;
};

Ray lunar_ray(const Context& c, ModelRef model, Vector observer_date = {}) {
    const Vector observer = c.frame.to_ecliptic(observer_date);
    const Vector receiver = c.earth_position_km + observer;
    Vector moon = c.moon_position_km;
    double delay = norm(moon - observer) / c_km_day;
    Vector ray = moon - observer;
    Vector heliocentric=c.earth_position_km+moon;
    auto previous_tdb = c.tdb;
    auto previous_tt = c.tt;
    // Retard BOTH the Earth and its geocentric Moon. Retarding only the Moon
    // then applying annual aberration would introduce a spurious ~20" term.
    for (int i = 0; i < 3; ++i) {
        const auto emission_tdb = shifted(c.tdb, -delay);
        const auto emission_tt = shifted(c.tt, -delay);
        // An exactly repeated pair of split instants would recompute the
        // same ray. No tolerance, truncated series or approximate cache.
        if (emission_tdb.day == previous_tdb.day && emission_tdb.frac == previous_tdb.frac
            && emission_tt.day == previous_tt.day && emission_tt.frac == previous_tt.frac)
            break;
        moon = model_position(model, emission_tdb);
        const Vector earth_emitted = earth_position(emission_tt);
        heliocentric=earth_emitted+moon;
        // Subtract the large Earth positions before adding the lunar vector
        // to retain lunar distance precision.
        ray = moon + (earth_emitted - receiver);
        delay = norm(ray) / c_km_day;
        previous_tdb = emission_tdb;
        previous_tt = emission_tt;
    }
    const Vector rotational_velocity_date = {
        -earth_rotation_rad_s * seconds_per_day * observer_date[1],
         earth_rotation_rad_s * seconds_per_day * observer_date[0], 0.0};
    const Vector velocity = c.earth_velocity_km_day
        + c.frame.to_ecliptic(rotational_velocity_date);
    const double distance = norm(ray);
    if (!(distance > moon_radius_km) || !std::isfinite(distance))
        throw std::invalid_argument("spp::moon: model/observer produces an invalid lunar range");
    const double travel=static_cast<double>(c.tdb.day-previous_tdb.day)+(c.tdb.frac-previous_tdb.frac);
    // Reuse the already evaluated reception/emission positions for the
    // short lunar velocity secant. No additional provider evaluations. The
    // Earth velocity is the same centered eleven-minute estimate used above.
    const Vector lunar_velocity=c.earth_velocity_km_day+(1.0/travel)*(c.moon_position_km-moon);
    return {aberration(ray, velocity),ray,heliocentric,lunar_velocity,velocity,distance,previous_tdb};
}

Vector illumination_source(const Ray& r) {
    // Illumination is determined at reflection, using the astrometric
    // observer-to-Moon ray, not the aberrated direction seen by the receiver.
    // The heliocentric approximation holds the Sun fixed; solar barycentric
    // motion during the incoming light time is outside the current model.
    const double solar_distance=norm(r.heliocentric_at_emission),lunar_speed=norm(r.moon_velocity_km_day);
    if(!(solar_distance>0) || !std::isfinite(solar_distance) || !std::isfinite(lunar_speed)
        || lunar_speed>=c_km_day)
        throw std::invalid_argument("spp::moon: invalid solar distance or lunar illumination velocity");
    return aberration((-1.0)*r.heliocentric_at_emission,r.moon_velocity_km_day);
}
double physical_phase(const Ray& r) {
    return angle((-1.0)*r.astrometric_j2000,illumination_source(r));
}

struct PhaseGeometry {
    Vector moon_ecliptic;
    Ecliptic solar;
    double cycle_deg;
};

PhaseGeometry phase_geometry(const Context& c, const Ray& r) {
    const Vector equ = c.frame.to_equatorial(r.apparent_j2000);
    const Vector ecl = rotate_x(equ, -c.nutation.true_obliquity_deg * deg2rad);
    const auto solar = sun::apparent_ecliptic(c.tt, c.nutation);
    const double lon = norm360(std::atan2(ecl[1], ecl[0]) * rad2deg);
    return {ecl, solar, norm360(lon - solar.lon_deg)};
}

double cycle(const Time& t, ModelRef model) {
    const Context c(t, model);
    return phase_geometry(c, lunar_ray(c, model)).cycle_deg;
}

Phase phase_from(const Context& c, const Ray& r) {
    const auto g = phase_geometry(c, r);
    const auto& ecl = g.moon_ecliptic;
    const auto& solar = g.solar;
    const Vector s = spherical(solar.lon_deg, solar.lat_deg, 1.0);
    const double elongation = angle(ecl,s) * rad2deg;
    // Keep the apparent elongation/cycle separate from the physical
    // illumination geometry. Apparent directions alone form a different
    // triangle because observer aberration changes its angles.
    const double i = physical_phase(r);
    const double half_cos = std::cos(0.5*i);
    return {g.cycle_deg, elongation, i * rad2deg, half_cos*half_cos, g.cycle_deg < 180.0};
}

struct ObserverSky {
    Equatorial equatorial;
    Horizontal horizontal;
    double geometric_alt_deg;
    double distance_km;
    double lst_deg;
    Ray ray;
};

ObserverSky observer_sky(const Context& c, const Time& t, const Observer& obs,
                         const Atmosphere& atm, ModelRef model) {
    const double lst = earth::lst_deg(t.ut1(), c.nutation, norm360(obs.lon_deg));
    const auto ray = lunar_ray(c, model, sky::observer_position_km(obs, lst));
    const auto eq = equatorial(c.frame.to_equatorial(ray.apparent_j2000));
    const auto hor = sky::to_horizontal(eq, obs.lat_deg, lst);
    return {eq, {sky::refraction(hor.alt_deg, sky::at_elevation(atm, obs.elev_m)), hor.az_deg},
            hor.alt_deg, ray.distance_km, lst,ray};
}

Disk disk_from(const Context& c, const Observer& obs, const ObserverSky& s,Vector incoming_sun) {
    // Apparent sky bearing (Meeus 48.5); illumination is a separate
    // reflection-geometry quantity and does not use this apparent triangle.
    const auto solar=sky::diurnal_aberration(
        sky::parallax(sun::apparent(c.tt,c.nutation),obs,s.lst_deg),obs,s.lst_deg);
    const double d=(solar.ra_deg-s.equatorial.ra_deg)*deg2rad;
    const double ds=solar.dec_deg*deg2rad,dm=s.equatorial.dec_deg*deg2rad;
    const double east=std::cos(ds)*std::sin(d);
    const double north=std::sin(ds)*std::cos(dm)-std::cos(ds)*std::sin(dm)*std::cos(d);
    const double i=angle((-1.0)*s.ray.astrometric_j2000,incoming_sun);
    const double half_cos=std::cos(0.5*i);
    const double phi=obs.lat_deg*deg2rad;
    const double h=wrap180(s.lst_deg-s.equatorial.ra_deg)*deg2rad;
    const double ze=std::cos(phi)*std::sin(h);
    const double zn=std::sin(phi)*std::cos(dm)-std::cos(phi)*std::sin(dm)*std::cos(h);
    const bool bright=std::hypot(east,north)>1e-12, zenith=std::hypot(ze,zn)>1e-12;
    return {i*rad2deg,half_cos*half_cos,bright?norm360(std::atan2(east,north)*rad2deg):0.0,
        zenith?wrap180(std::atan2(ze,zn)*rad2deg):0.0,bright,zenith};
}

Surface surface_from(const Context& c,const Ray& ray,const Equatorial& eq,RotationRef spin,Vector incoming_sun) {
    const auto rotation=spin(ray.emission_tdb);
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
        if(!std::isfinite(rotation.matrix[i][j]) || !std::isfinite(rotation.derivative_per_day[i][j])
            || std::abs(dot(rotation.matrix[i],rotation.matrix[j])-(i==j?1:0))>1e-8)
            throw std::invalid_argument("spp::moon: rotation provider must return a proper orthonormal frame");
    }
    if(dot(cross(rotation.matrix[0],rotation.matrix[1]),rotation.matrix[2])<0)
        throw std::invalid_argument("spp::moon: rotation provider returned a reflected frame");
    const auto to_body=rotation*spp::frames::mean_ecliptic_j2000_to_icrf();
    using Direction=CartesianPosition<origin::Selenocentric,frame::MeanEclipticJ2000>;
    const auto observer=to_body(Direction{(-1.0)*ray.astrometric_j2000}).position_km;
    const auto solar=to_body(Direction{incoming_sun}).position_km;
    auto longitude=[](Vector v){return wrap180(std::atan2(v[1],v[0])*rad2deg);};
    auto latitude=[](Vector v){return std::atan2(v[2],std::hypot(v[0],v[1]))*rad2deg;};
    Surface result{longitude(observer),latitude(observer),longitude(solar),latitude(solar),0,false,{}};
    const double r=eq.ra_deg*deg2rad,d=eq.dec_deg*deg2rad;
    const Vector east{-std::sin(r),std::cos(r),0};
    const Vector north{-std::sin(d)*std::cos(r),-std::sin(d)*std::sin(r),std::cos(d)};
    const Vector center{std::cos(d)*std::cos(r),std::cos(d)*std::sin(r),std::sin(d)};
    const Vector n=(1.0/ray.distance_km)*ray.astrometric_j2000;
    const Vector apparent=(1.0/ray.distance_km)*ray.apparent_j2000;
    const Vector beta=(1.0/c_km_day)*ray.observer_velocity_km_day;
    const double scale=std::sqrt(1-dot(beta,beta))/(1+dot(n,beta));
    for(std::size_t j=0;j<3;++j) {
        // The aberration map is conformal on the sky. Divide out its common
        // angular scale to obtain a proper orientation of an orthographic
        // disk, preserving its physical sub-observer point and pole flags.
        const Vector axis0=to_body.matrix[j];
        const double radial=dot(n,axis0);
        const Vector tangent=axis0-radial*n;
        const Vector mapped=(1.0/scale)*aberration_tangent(n,apparent,tangent,ray.observer_velocity_km_day)
            +radial*apparent;
        const auto axis=c.frame.to_equatorial(mapped);
        result.body_to_view[0][j]=-dot(east,axis);
        result.body_to_view[1][j]=dot(north,axis);
        result.body_to_view[2][j]=-dot(center,axis);
    }
    result.pole_defined=std::hypot(result.body_to_view[0][2],result.body_to_view[1][2])>1e-12;
    if(result.pole_defined)result.pole_position_angle_deg=norm360(std::atan2(
        -result.body_to_view[0][2],result.body_to_view[1][2])*rad2deg);
    return result;
}

}  // namespace

Observation observe(const Time& t,const Observer& obs,const Atmosphere& atm,ModelRef model,RotationRef spin) {
    validate(obs, atm);
    const Context c(t, model);
    const auto geocentric = lunar_ray(c, model);
    const auto s = observer_sky(c, t, obs, atm, model);
    const auto incoming_sun=illumination_source(s.ray);
    return {s.equatorial, s.horizontal, s.geometric_alt_deg, s.distance_km,
            std::asin(moon_radius_km/s.distance_km)*rad2deg,
            phase_from(c, geocentric), disk_from(c,obs,s,incoming_sun),surface_from(c,s.ray,s.equatorial,spin,incoming_sun)};
}
Observation observe(const Time& t,const Observer& obs,const Atmosphere& atm,ModelRef model) {
    const spp::model::iau_moon::Rotation spin;
    return observe(t,obs,atm,model,RotationRef(spin));
}
Surface surface(const Time& t,const std::optional<Observer>& obs,ModelRef model,RotationRef spin) {
    if(obs)validate(*obs,{});
    const Context c(t,model);
    if(obs) {
        const auto s=observer_sky(c,t,*obs,{.refraction=false},model);
        return surface_from(c,s.ray,s.equatorial,spin,illumination_source(s.ray));
    }
    const auto ray=lunar_ray(c,model);
    const auto eq=equatorial(c.frame.to_equatorial(ray.apparent_j2000));
    return surface_from(c,ray,eq,spin,illumination_source(ray));
}
Surface surface(const Time& t,const std::optional<Observer>& obs,ModelRef model) {
    const spp::model::iau_moon::Rotation spin;
    return surface(t,obs,model,RotationRef(spin));
}

Horizontal horizontal(const Time& t, const Observer& obs, const Atmosphere& atm, ModelRef model) {
    validate(obs,atm);
    const Context c(t,model);
    return observer_sky(c,t,obs,atm,model).horizontal;
}

Coordinates coordinates(const Time& t,const Observer& obs,const Atmosphere& atm,ModelRef model) {
    validate(obs, atm);
    const Context c(t, model);
    const auto s=observer_sky(c,t,obs,atm,model);
    return {s.equatorial,s.horizontal,s.geometric_alt_deg,s.distance_km,
        std::asin(radius_km/s.distance_km)*rad2deg,wrap180(s.lst_deg-s.equatorial.ra_deg)};
}

Equatorial apparent(const Time& t, ModelRef model) {
    const Context c(t, model);
    return equatorial(c.frame.to_equatorial(lunar_ray(c, model).apparent_j2000));
}

Phase phase(const Time& t, ModelRef model) {
    const Context c(t, model);
    return phase_from(c, lunar_ray(c, model));
}

std::vector<PhaseEvent> phases(const Time& start, const Time& end, double tolerance_s, ModelRef model) {
    if (!std::isfinite(tolerance_s) || tolerance_s < 1.0e-5 || tolerance_s > 60.0)
        throw std::invalid_argument("spp::moon::phases: tolerance must be between 0.00001 and 60 seconds");
    const auto a = start.utc(), b = end.utc();
    if (b.day < a.day || (b.day == a.day && b.frac < a.frac))
        throw std::invalid_argument("spp::moon::phases: reversed interval");
    if (a.day == b.day && a.frac == b.frac) return {};
    // Reject an unsupported endpoint before doing a potentially long search.
    (void)model(start.tdb());
    (void)model(end.tdb());
    // One-day brackets are specific to the terrestrial Moon's monotonic
    // synodic longitude. They are not a general satellite-event algorithm.
    Time current = start;
    double cycle_lo = cycle(current, model);
    double target = 90.0 * std::ceil(cycle_lo / 90.0);
    std::vector<PhaseEvent> result;
    while (true) {
        const auto u = current.utc();
        // Subtract integer days only when they differ by at most one. This
        // keeps a distant end date out of both integer overflow and the local
        // fractional-day arithmetic. Every refinement stays within one day.
        const double remaining = b.day > u.day + 1 ? 2.0
            : static_cast<double>(b.day-u.day) + (b.frac-u.frac);
        const double width = std::min(1.0, remaining);
        const Time next = current.shifted(width);
        const double wrapped_hi = cycle(next, model);
        const double cycle_hi = cycle_lo + wrap180(wrapped_hi - norm360(cycle_lo));
        if (!(cycle_hi > cycle_lo) || cycle_hi-cycle_lo > 45.0)
            throw std::runtime_error("spp::moon::phases: model violates terrestrial Moon search assumptions");
        if (target >= cycle_lo && target <= cycle_hi) {
            auto f = [&](double offset) {
                return wrap180(cycle(current.shifted(offset), model) - norm360(target));
            };
            const double root = target == cycle_lo ? 0.0 : (target == cycle_hi ? width
                : refine_root(f, width, cycle_hi-target, 0.0, cycle_lo-target, tolerance_s/seconds_per_day, 80));
            if (std::abs(f(root)) > 20.0 * tolerance_s/seconds_per_day + 1.0e-8)
                throw std::runtime_error("spp::moon::phases: root failed to converge");
            const Time time = current.shifted(root);
            const auto instant = time.utc();
            if (instant.day < b.day || (instant.day == b.day && instant.frac < b.frac)) {
                const int index = static_cast<int>(std::lround(std::fmod(target, 360.0)/90.0)) % 4;
                result.push_back({static_cast<PrincipalPhase>(index), time});
            }
            target += 90.0;
        }
        if (remaining <= 1.0) break;
        current = next;
        cycle_lo = cycle_hi;
    }
    return result;
}

}  // namespace spp::moon::detail
