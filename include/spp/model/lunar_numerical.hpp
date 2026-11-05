#pragma once

#include "spp/model/integrator.hpp"
#include "spp/model/spk.hpp"
#include "spp/core/rotation.hpp"

namespace spp::model::lunar_numerical {
using State=CartesianState<origin::Geocentric,frame::Icrf>;
using Vector=std::array<double,3>;
struct EarthTideParameters {
    std::array<double,3> love{0.335,0.277,0.283};
    std::array<double,3> orbit_lag_days{0.077,-0.01,-0.20982243826832672};
    std::array<double,3> rotation_lag_days{0,0.011929122823877609,-0.00074236275298830701};
};
struct ForceOptions {
    bool relativity=true;       // major-body PPN, beta=gamma=1
    bool earth_oblateness=true; // zonal field with pole of date
    bool lunar_figure=true;     // PA field, IAU orientation + DE440 ME/PA offset
    double earth_j2=1.08262539e-3;
    double earth_radius_km=6378.1366;
    double moon_j2=2.03214360135e-4;
    double moon_c22=2.2380845245752079e-5;
    double moon_radius_km=1738;
    unsigned earth_degree=5;   // zonal field, degree 2..5
    unsigned moon_degree=6;    // fixed PA field, degree/order 2..6
    bool earth_tides=true;     // Sun- and Moon-raised, order-dependent response
    EarthTideParameters earth_tide;
    bool lunar_tides=true;     // delayed Earth-raised tide and spin deformation
    double moon_love=0.02419;
    double moon_lag_days=0.14519388258636395;
    double moon_reference_spin_rad_day=13.17635815*deg2rad;
    bool external_figure_interactions=true; // figure/point-mass reactions through Saturn
    bool solar_oblateness=true;
    bool solar_radiation_pressure=true;
    bool solar_frame_dragging=true; // also requires relativity=true
    double earth_srp_fraction=2e-14;
    double moon_srp_fraction=1.44e-13;
    // Cache smooth, time-only external states/orientations during propagation.
    // Disable for direct-reference comparisons. Never interpolates the Moon's
    // propagated state and never changes the force terms or RK tolerances.
    bool environment_cache=true;
};

// Geocentric translational model: Earth/Moon attraction, prescribed DE440 EMB,
// Sun and planetary-system perturbations, optional major-body 1PN, Earth
// zonals, lunar harmonics, delayed tidal/spin deformation and small solar
// figure/radiation/frame-dragging terms. Does not include
// asteroids or integrated rotational dynamics;
// it is not a reimplementation of DE440. Configuration is explicit.
class Forces {
public:
    explicit Forces(spk::Kernel kernel,ForceOptions options={});
    Vector operator()(Instant<TDB> t,const State& state)const;
    const spk::Kernel& kernel()const{return kernel_;}
    Forces prepared(Instant<TDB> first,Instant<TDB> last)const;
    std::size_t environment_bytes()const;
    std::size_t environment_samples()const;
private:
    struct Environment;
    spk::Kernel kernel_;ForceOptions options_;
    std::shared_ptr<const Environment> environment_;
};

class Ephemeris {
public:
    // Seed from this kernel's geometric Moon/Earth state once. Forces use
    // planets/EMB, never the kernel's Earth/Moon states after construction.
    Ephemeris(spk::Kernel kernel,Instant<TDB> seed,Instant<TDB> first,Instant<TDB> last,
              IntegrationOptions integration={},ForceOptions forces={});
    // Explicit ICRF/geocentric seed, for a fitted state or alternate ephemeris.
    Ephemeris(Forces forces,Instant<TDB> seed,State state,Instant<TDB> first,Instant<TDB> last,
              IntegrationOptions integration={});
    GeocentricEclipticJ2000State operator()(Instant<TDB> t)const;
    GeocentricEclipticJ2000Position position(Instant<TDB> t)const;
    State icrf(Instant<TDB> t)const{return orbit_(t);}
    const IntegrationStats& stats()const{return orbit_.stats();}
    std::size_t storage_bytes()const{return orbit_.storage_bytes();}
    std::size_t environment_bytes()const{return orbit_.force().environment_bytes();}
    std::size_t environment_samples()const{return orbit_.force().environment_samples();}
private:
    IntegratedOrbit<Forces> orbit_;
};
} // namespace spp::model::lunar_numerical
