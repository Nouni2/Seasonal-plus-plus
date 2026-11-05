#include "spp/model/lunar_numerical.hpp"
#include "spp/moon.hpp"
#include "spp/core/angles.hpp"
#include "core/vector.hpp"
#include "core/gravity.hpp"
#include "core/tides.hpp"
#include "core/chebyshev.hpp"
#include "model/lunar_gravity.hpp"
#include "phase_reference.hpp"
#include "error_stats.hpp"
#include "model/solar_effects.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdexcept>

namespace {
using namespace spp::detail;
using Epoch=spp::Instant<spp::TDB>;
using State=spp::model::lunar_numerical::State;
using Clock=std::chrono::steady_clock;
void require(bool condition,const char* why){if(!condition)throw std::runtime_error(why);}
Epoch shift(Epoch t,double d){return {t.day,t.frac+d};}
struct Central {
    static constexpr double mu=403503.2356254802*spp::seconds_per_day*spp::seconds_per_day;
    Vector operator()(Epoch,const State& state)const {
        const auto& r=state.position_km;return (-mu/std::pow(norm(r),3))*r;
    }
};
void gravity_controls() {
    const auto polynomial=ChebyshevSeries<3,12>::fit([](double x){return Vector{1+2*x-3*x*x,1e8+300*x*x*x,0.1*std::pow(x,12)};});
    for(int i=0;i<=1024;++i) {
        const double x=-1+2.0*i/1024;const auto v=polynomial(x);
        require(std::abs(v[0]-(1+2*x-3*x*x))<1e-13 && std::abs(v[1]-(1e8+300*x*x*x))<6e-8
            && std::abs(v[2]-0.1*std::pow(x,12))<1e-15,"Chebyshev interpolation fails exact-polynomial/endpoints control");
    }
    for(double delta:{-0.9999999999999998,-0.9,-0.01,-1e-10,-1e-25,0.0,1e-25,1e-10,0.01,0.9,1e6}) {
        const double fast=delta*(3+delta*(3+delta))/(1+(1+delta)*std::sqrt(1+delta));
        const long double exact=std::expm1(1.5L*std::log1p(static_cast<long double>(delta)));
        require(std::abs(static_cast<long double>(fast)-exact)<8*std::numeric_limits<double>::epsilon()*std::max(1e-30L,std::abs(exact)),
            "rationalized differential gravity loses precision");
    }
    std::ifstream file(SPP_GRAVITY_REFERENCE);require(file.good(),"cannot open independent gravity fixture");
    std::size_t count=0;double maximum=0,relative=0,zonal_error=0,checksum=0;
    for(std::string line;std::getline(file,line);) {
        if(line.empty() || line.front()=='#')continue;
        for(auto& c:line)if(c==',')c=' ';
        std::istringstream in(line);std::string body;unsigned degree;Vector r,reference;double radius,mu;
        in>>body>>degree>>r[0]>>r[1]>>r[2]>>radius>>mu>>reference[0]>>reference[1]>>reference[2];
        require(!in.fail(),"invalid gravity fixture");
        auto c=spp::model::lunar_numerical::detail::moon_cosine,s=spp::model::lunar_numerical::detail::moon_sine;
        if(body=="earth") {
            c={};s={};c[2][0]=-1.08262539e-3;c[3][0]=2.53241e-6;c[4][0]=1.619898e-6;c[5][0]=2.277345e-7;
        } else require(body=="moon","unknown gravity reference body");
        const auto actual=harmonic_acceleration<6>(r,radius,mu,degree,c,s);
        if(degree==2) {
            const auto fast=quadrupole_acceleration(r,radius,mu,c[2][0],c[2][1],s[2][1],c[2][2],s[2][2]);
            require(norm(fast-reference)<1e-20+64*std::numeric_limits<double>::epsilon()*norm(reference),"quadrupole reaction field differs from oracle");
        }
        const double error=norm(actual-reference);maximum=std::max(maximum,error);
        relative=std::max(relative,error/std::max(1e-30,norm(reference)));
        require(error<1e-20+64*std::numeric_limits<double>::epsilon()*norm(reference),"harmonic potential-gradient oracle differs");
        if(body=="earth") {
            const std::array<double,6> j{0,0,1.08262539e-3,-2.53241e-6,-1.619898e-6,-2.277345e-7};
            const auto fast=zonal_acceleration(r,{0,0,1},radius,mu,degree,j);
            zonal_error=std::max(zonal_error,norm(fast-actual));
            require(norm(fast-reference)<1e-20+64*std::numeric_limits<double>::epsilon()*norm(reference),"fast zonal field differs from oracle");
        }
        checksum+=actual[0]+actual[1]+actual[2];++count;
    }
    require(count==28,"gravity reference coverage changed");
    const Vector sun_r{1.4e8,-3e7,2e7},sun_v{1e6,2e6,-3e5};
    const auto pressure=spp::model::de440_solar::radiation(sun_r,Central::mu,1e-13);
    require(dot(pressure,sun_r)>0 && norm(cross(pressure,sun_r))<1e-20,"radiation pressure does not point away from the Sun");
    const auto magnetic=spp::model::de440_solar::frame_dragging(sun_r,sun_v,Central::mu);
    require(std::abs(dot(magnetic,sun_v))<1e-20,"solar frame dragging does work along velocity");
    require(norm(magnetic+spp::model::de440_solar::frame_dragging(sun_r,(-1.0)*sun_v,Central::mu))<1e-30,
        "frame dragging velocity symmetry failed");
    namespace solar=spp::model::de440_solar;
    const auto axis=solar::pole();const double distance=1.5e8,gm_sun=132712440041.27942*std::pow(spp::seconds_per_day,2);
    const Vector velocity{2e6,0,0};
    const double spin_magnitude=gm_sun*solar::radius_km*solar::radius_km*solar::inertia_fraction*solar::spin_rad_day;
    const auto expected_lt=(4*spin_magnitude/(std::pow(299792.458*spp::seconds_per_day,2)*std::pow(distance,3)))*cross(axis,velocity);
    require(norm(solar::frame_dragging(distance*axis,velocity,gm_sun)-expected_lt)<1e-21,"polar Lense-Thirring magnitude differs");
    const auto expected_j2=(3*gm_sun*solar::j2*solar::radius_km*solar::radius_km/std::pow(distance,4))*axis;
    require(norm(solar::oblateness(distance*axis,gm_sun)-expected_j2)<1e-20,"solar J2 polar acceleration differs");
    std::ifstream tides(std::filesystem::path(SPP_GRAVITY_REFERENCE).parent_path()/"earth_tides.csv");
    require(tides.good(),"cannot open independent tidal fixture");
    std::size_t tidal_count=0;double tidal_oracle=0;
    for(std::string line;std::getline(tides,line);) {
        if(line.empty() || line.front()=='#')continue;
        for(char& c:line)if(c==',')c=' ';
        std::istringstream in(line);unsigned order;Vector r,q,expected;double love,radius,mu;
        in>>order>>r[0]>>r[1]>>r[2]>>q[0]>>q[1]>>q[2]>>love>>radius>>mu>>expected[0]>>expected[1]>>expected[2];
        require(!in.fail(),"invalid tidal fixture");
        const auto actual=tide_acceleration(r,q,{0,0,1},order,radius,mu,love);
        tidal_oracle=std::max(tidal_oracle,norm(actual-expected));
        require(norm(actual-expected)<1e-27+64*std::numeric_limits<double>::epsilon()*norm(expected),"individual tidal mode differs from Decimal potential gradient");
        ++tidal_count;
    }
    require(tidal_count==18,"tidal oracle coverage changed");
    double tide_error=0;
    for(Vector pole:{Vector{0,0,1},Vector{0.6,0,0.8}})for(Vector r:{Vector{384400,0,0},Vector{1e5,-3e5,2e5}}) {
        const Vector q{2e5,-1e5,3e5};const double radius=6378.1366,mu=4902.80011845755*std::pow(spp::seconds_per_day,2),love=0.3;
        Vector actual{};for(unsigned order=0;order<3;++order)actual=actual+tide_acceleration(r,q,pole,order,radius,mu,love);
        const double r2=dot(r,r),q2=dot(q,q),rq=dot(r,q);
        const auto expected=(1.5*love*mu*std::pow(radius,5)/(std::pow(norm(r),5)*std::pow(norm(q),5)))
            *(2*rq*q+(q2-5*rq*rq/r2)*r);
        const double e=norm(actual-expected);tide_error=std::max(tide_error,e);
        require(e<1e-21,"tidal modes do not sum to the independent equilibrium field");
    }
    std::printf("Gravity gradients: 28 independent Decimal potential derivatives; max %.9g km/day2, relative %.9g, zonal %.9g; tidal tensor control %.9g km/day2; checksum %.9g\n",
        maximum,relative,zonal_error,tide_error,checksum);
    std::printf("Tidal modes: 18 independent potential gradients, max %.9g km/day2\n",tidal_oracle);
    const auto begin=Clock::now();
    for(int i=0;i<16384;++i) {
        const double a=i*0.001;const Vector r{384400*std::cos(a),384400*std::sin(a),7000};
        const auto g=harmonic_acceleration<6>(r,1738,Central::mu,6,
            spp::model::lunar_numerical::detail::moon_cosine,spp::model::lunar_numerical::detail::moon_sine);
        checksum+=g[0]+g[1]+g[2];
    }
    std::printf("Degree/order 6 field: %.3f us/query, checksum %.9g\n",
        std::chrono::duration<double,std::micro>(Clock::now()-begin).count()/16384,checksum);
    require(std::isfinite(checksum),"non-finite gravity benchmark checksum");
}
void numerical() {
    constexpr double radius=384400,omega=std::sqrt(Central::mu/(radius*radius*radius));
    const Epoch seed{2451545,0.3};const State state{{radius,0,0},{0,radius*omega,0}};
    const auto start=Clock::now();
    const spp::model::IntegratedOrbit<Central> orbit(Central{},seed,state,shift(seed,-32),shift(seed,32));
    double maxp=0,maxv=0,derivative=0,energy=0,checksum=0;
    for(int i=0;i<=4096;++i) {
        const double d=-32+64.0*i/4096;
        const auto s=orbit(shift(seed,d));
        const Vector p{radius*std::cos(omega*d),radius*std::sin(omega*d),0};
        const Vector v{-radius*omega*std::sin(omega*d),radius*omega*std::cos(omega*d),0};
        maxp=std::max(maxp,norm(s.position_km-p));maxv=std::max(maxv,norm(s.velocity_km_per_day-v));
        energy=std::max(energy,std::abs((0.5*dot(s.velocity_km_per_day,s.velocity_km_per_day)-Central::mu/norm(s.position_km))
            /(-Central::mu/(2*radius))-1));
        if(i>0 && i<4096) {
            constexpr double h=1e-4;
            const auto difference=(1.0/(12*h))*((-1.0)*orbit.position(shift(seed,d+2*h)).position_km
                +8*orbit.position(shift(seed,d+h)).position_km-8*orbit.position(shift(seed,d-h)).position_km
                +orbit.position(shift(seed,d-2*h)).position_km);
            derivative=std::max(derivative,norm(difference-s.velocity_km_per_day));
        }
    }
    const double elapsed=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    std::printf("Circular orbit/dense output +/-32 days: %.9g m, %.9g km/day, derivative %.9g km/day, energy %.9g; %.3f ms incl controls\n",
        maxp*1000,maxv,derivative,energy,elapsed);
    require(maxp<1e-5 && maxv<1e-5 && derivative<1e-4 && energy<1e-10,"integrator analytic accuracy gate failed");
    const auto a=orbit.stats();
    require(a.accepted_steps>0 && a.force_evaluations==2+6*(a.accepted_steps+a.rejected_steps),"integrator FSAL/counter contract failed");
    const auto begin=Clock::now();
    for(int repeat=0;repeat<64;++repeat)for(int i=0;i<4096;++i) {
        const auto s=orbit(shift(seed,-32+64.0*((i*1543)%4096)/4096));checksum+=s.position_km[0]+s.velocity_km_per_day[1];
    }
    const double us=std::chrono::duration<double,std::micro>(Clock::now()-begin).count()/(64*4096);
    std::printf("Dense random query %.3f us; %zu accepted, %zu rejected, %zu force calls, %zu bytes; checksum %.12g\n",
        us,a.accepted_steps,a.rejected_steps,a.force_evaluations,orbit.storage_bytes(),checksum);
    require(std::isfinite(checksum) && orbit.stats().force_evaluations==a.force_evaluations,"dense queries evaluated force");
    const auto end=orbit(shift(seed,32));
    const spp::model::IntegratedOrbit<Central> reversed(Central{},shift(seed,32),end,seed,shift(seed,32));
    require(norm(reversed(seed).position_km-state.position_km)<2e-5,"backward propagation fails roundtrip");
    bool out=false;try{(void)orbit(shift(seed,32.001));}catch(const std::out_of_range&){out=true;}
    require(out,"integrator extrapolates outside explicit window");
    for(int bad=0;bad<4;++bad) {
        spp::model::IntegrationOptions options;
        if(bad==0)options.relative_tolerance=std::numeric_limits<double>::quiet_NaN();
        if(bad==1)options.minimum_step_days=1;
        if(bad==2)options.maximum_steps=0;
        bool rejected=false;
        try {const spp::model::IntegratedOrbit<Central> invalid(Central{},seed,state,
            bad==3?shift(seed,1):shift(seed,-1),shift(seed,1),options);}
        catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"invalid integration configuration accepted");
    }
    spp::model::IntegrationOptions adaptive;adaptive.maximum_step_days=1;adaptive.initial_step_days=1;
    const spp::model::IntegratedOrbit<Central> automatic(Central{},seed,state,seed,shift(seed,1),adaptive);
    require(automatic.stats().rejected_steps>0,"adaptive benchmark did not exercise error rejection");
    require(norm(automatic(shift(seed,1)).position_km-orbit(shift(seed,1)).position_km)<1e-5,
        "adaptive rejection changed solution beyond numerical budget");
    adaptive.maximum_steps=1;
    bool exhausted=false;
    try{const spp::model::IntegratedOrbit<Central> limited(Central{},seed,state,seed,shift(seed,1),adaptive);}
    catch(const std::runtime_error&){exhausted=true;}
    require(exhausted,"step budget ignored");
}
void environment_controls(const spp::model::spk::Kernel& kernel) {
    using spp::model::lunar_numerical::Forces;
    const Forces direct(kernel);double error=0,relative=0;std::size_t count=0;
    for(Epoch seed:{Epoch{2415020,0.5},Epoch{2433282,0.5},Epoch{2451545,0},Epoch{2460409,0.25},
                    Epoch{2461348,0},Epoch{2469807,0.5},Epoch{2488069,0.5}}) {
        const auto first=shift(seed,-6.17),last=shift(seed,6.13);const auto cached=direct.prepared(first,last);
        require(cached.environment_bytes()>0 && direct.environment_bytes()==0,"environment ownership/cache controls failed");
        // Endpoints, partial final segments, boundaries and nonuniform trial
        // states: the cache must not depend on the propagated state.
        for(int i=0;i<=512;++i) {
            const auto t=shift(first,12.3*((i*197)%513)/512.0);
            const auto s=kernel.state(301,399,t);
            const State state{s.position_km+Vector{0.13*i,-0.07*i,0.03*i},s.velocity_km_per_day};
            const auto a=direct(t,state),b=cached(t,state);
            error=std::max(error,norm(a-b));relative=std::max(relative,norm(a-b)/norm(a));++count;
        }
        for(const auto t:{first,shift(first,4),shift(first,8),last}) {
            const auto s=kernel.state(301,399,t);const State state{s.position_km,s.velocity_km_per_day};
            error=std::max(error,norm(direct(t,state)-cached(t,state)));++count;
        }
        for(const auto t:{shift(first,-0.1),shift(last,0.1)}) {
            const auto s=kernel.state(301,399,t);const State state{s.position_km,s.velocity_km_per_day};
            require(direct(t,state)==cached(t,state),"environment cache extrapolated outside its window");
        }
    }
    bool reversed=false;try{(void)direct.prepared({2451546,0},{2451545,0});}catch(const std::invalid_argument&){reversed=true;}
    require(reversed,"reversed environment domain accepted");
    bool coverage=false;
    try{(void)direct.prepared({2396752,0},{2396753,0});}catch(const std::out_of_range&){coverage=true;}
    require(coverage,"cached environment hid missing kernel endpoint coverage");
    spp::model::lunar_numerical::ForceOptions options;options.environment_cache=false;
    require(Forces(kernel,options).prepared({2451545,0},{2451546,0}).environment_bytes()==0,"direct environment option ignored");
    std::printf("Environment interpolation: %zu force comparisons at seven epochs, max %.9g km/day2, relative %.9g; domain/endpoints/direct mode PASS\n",count,error,relative);
    require(error<1e-10,"environment interpolation exceeds sampled force budget");
}
void physical(const char* path) {
    using spp::model::lunar_numerical::Ephemeris;
    static_assert(spp::moon::Ephemeris<Ephemeris>);
    const spp::model::spk::Kernel kernel(path);
    environment_controls(kernel);
    double maximum=0,checksum=0;
    for(Epoch seed:{Epoch{2415020,0.5},Epoch{2433282,0.5},Epoch{2451545,0},Epoch{2460409,0.25},
                    Epoch{2461348,0},Epoch{2469807,0.5},Epoch{2488069,0.5}}) {
        const auto begin=Clock::now();const Ephemeris orbit(kernel,seed,shift(seed,-16),shift(seed,16));
        const double ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
        double maxp=0,maxv=0;
        for(int i=0;i<=256;++i) {
            const auto t=shift(seed,-16+i/8.0);const auto a=orbit.icrf(t);const auto b=kernel.state(301,399,t);
            maxp=std::max(maxp,norm(a.position_km-b.position_km));
            maxv=std::max(maxv,norm(a.velocity_km_per_day-b.velocity_km_per_day));
            checksum+=a.position_km[0];
        }
        maximum=std::max(maximum,maxp);
        const auto stats=orbit.stats();
        std::printf("Major-body+1PN+Earth J2+lunar figure vs DE440s: seed JD %lld + %.3f, +/-16 d; %.6f km, %.6f km/day; build %.3f ms, %zu steps/%zu rejected/%zu force calls\n",
            static_cast<long long>(seed.day),seed.frac,maxp,maxv,ms,stats.accepted_steps,stats.rejected_steps,stats.force_evaluations);
        require(maxp<0.005 && maxv<0.002,"reduced lunar force model outside 16-day physical budget");
        spp::model::lunar_numerical::ForceOptions direct_options;direct_options.environment_cache=false;
        const auto direct_start=Clock::now();
        const Ephemeris direct(kernel,seed,shift(seed,-16),shift(seed,16),{},direct_options);
        const double direct_ms=std::chrono::duration<double,std::milli>(Clock::now()-direct_start).count();
        double cache_error=0;
        for(int i=0;i<=1024;++i) {
            const auto t=shift(seed,-16+32.0*i/1024);
            cache_error=std::max(cache_error,norm(direct.icrf(t).position_km-orbit.icrf(t).position_km));
        }
        std::printf("  Environment cached/direct: %.9g m trajectory difference, %.3f/%.3f ms build, %zu samples, %zu auxiliary bytes\n",
            cache_error*1000,ms,direct_ms,orbit.environment_samples(),orbit.environment_bytes());
        require(cache_error<1e-7,"environment cache changed short-arc trajectory by over 0.1 mm");
        if(seed.day==2461348) {
            spp::model::IntegrationOptions options;options.maximum_step_days=1.0/96;
            options.relative_tolerance=1e-13;options.position_tolerance_km=1e-8;options.velocity_tolerance_km_day=1e-8;
            const Ephemeris fine(kernel,seed,shift(seed,-16),shift(seed,16),options);
            double refinement=0;
            for(int i=0;i<=512;++i) {
                const auto t=shift(seed,-16+i/16.0);
                refinement=std::max(refinement,norm(fine.icrf(t).position_km-orbit.icrf(t).position_km));
            }
            std::printf("  Half maximum step + tighter tolerances: max position change %.9g m\n",refinement*1000);
            require(refinement<1e-6,"lunar integration refinement error exceeds physical comparison budget");
            // End-to-end observer path against the same kernel's unpropagated
            // geometric states. This isolates propagation from Spp's sky model.
            struct KernelMoon {
                const spp::model::spk::Kernel& kernel;
                spp::moon::Position position(Epoch t)const {
                    const auto s=kernel.state(301,399,t);
                    const spp::CartesianPosition<spp::origin::Geocentric,spp::frame::Icrf> p{s.position_km};
                    return spp::frames::mean_ecliptic_j2000_to_icrf().inverse()(p);
                }
            } reference{kernel};
            const spp::Observer paris{48.8566,2.3522,35};
            double sky_error=0,surface_error=0;
            std::vector<spp::Time> dates;
            for(int i=0;i<128;++i) {
                const auto time=spp::Time::from_utc(2026,11,3,12).shifted(-15.5+31.0*i/127);dates.push_back(time);
                const auto a=spp::moon::observe(time,paris,orbit),b=spp::moon::observe(time,paris,reference);
                sky_error=std::max(sky_error,std::hypot(a.geometric_alt_deg-b.geometric_alt_deg,
                    spp::wrap180(a.horizontal.az_deg-b.horizontal.az_deg))*3600);
                surface_error=std::max(surface_error,std::hypot(spp::wrap180(a.surface.subobserver_lon_deg-b.surface.subobserver_lon_deg),
                    a.surface.subobserver_lat_deg-b.surface.subobserver_lat_deg)*3600);
            }
            std::printf("  Propagated vs kernel observer directions: max %.9g arcsec, surface %.9g arcsec\n",sky_error,surface_error);
            require(sky_error<0.01 && surface_error<0.01,"propagated observation/reference disagreement");
            const auto force_calls=orbit.stats().force_evaluations;
            for(bool complete:{false,true}) {
                const auto started=Clock::now();double observed_checksum=0;
                for(int batch=0;batch<16;++batch)for(const auto& time:dates) {
                    if(complete) {
                        const auto q=spp::moon::observe(time,paris,orbit);
                        observed_checksum+=q.horizontal.alt_deg+q.horizontal.az_deg+q.surface.subobserver_lon_deg;
                    } else {
                        const auto q=spp::moon::coordinates(time,paris,orbit);
                        observed_checksum+=q.horizontal.alt_deg+q.horizontal.az_deg;
                    }
                }
                const double us=std::chrono::duration<double,std::micro>(Clock::now()-started).count()/(16*dates.size());
                std::printf("  Integrated %s: %.3f us/query resident, checksum %.12g\n",complete?"full observation":"coordinates",us,observed_checksum);
                require(std::isfinite(observed_checksum),"invalid integrated observation timing checksum");
            }
            require(orbit.stats().force_evaluations==force_calls,"observation queries reintegrated trajectory");
        }
    }
    require(std::isfinite(checksum),"invalid physical integrator checksum");
    std::printf("Reduced-model maximum %.6f km; checksum %.12g. This is not DE440's complete dynamical model.\n",maximum,checksum);
    const Epoch seed{2461348,0};
    const auto begin=Clock::now();const Ephemeris extended(kernel,seed,shift(seed,-90),shift(seed,90));
    const double construction_ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
    double long_error=0;
    for(int i=0;i<=1440;++i) {
        const auto t=shift(seed,-90+i/8.0);const auto a=extended.icrf(t);const auto b=kernel.state(301,399,t);
        long_error=std::max(long_error,norm(a.position_km-b.position_km));
    }
    std::printf("GUI +/-90-day window: max %.9g m vs DE440s, %.3f ms construction, %zu bytes trajectory\n",
        long_error*1000,construction_ms,extended.storage_bytes());
    require(long_error<0.02,"reduced lunar model outside 90-day sampled budget");
    spp::model::lunar_numerical::ForceOptions direct_options;direct_options.environment_cache=false;
    const auto direct_begin=Clock::now();
    const Ephemeris direct(kernel,seed,shift(seed,-90),shift(seed,90),{},direct_options);
    const double direct_ms=std::chrono::duration<double,std::milli>(Clock::now()-direct_begin).count();
    double cache_error=0;
    for(int i=0;i<=2880;++i) {
        const auto t=shift(seed,-90+i/16.0);
        cache_error=std::max(cache_error,norm(direct.icrf(t).position_km-extended.icrf(t).position_km));
    }
    std::printf("  +/-90-day environment cached/direct: %.9g m trajectory difference, %.3f/%.3f ms build, %zu samples, %zu auxiliary bytes\n",
        cache_error*1000,construction_ms,direct_ms,extended.environment_samples(),extended.environment_bytes());
    require(cache_error<1e-7,"environment cache changed 90-day trajectory by over 0.1 mm");
}
void ablation(const char* path) {
    using namespace spp::model::lunar_numerical;
    const spp::model::spk::Kernel kernel(path);const Epoch seed{2461348,0};
    struct Case {const char* name;ForceOptions options;};
    std::array<Case,8> cases{{{"quadrupole",{}},{"higher gravity",{}},{"+ Earth tides",{}},{"+ lunar tides",{}},{"+ figure reactions",{}},
        {"+ solar J2",{}},{"+ photon pressure",{}},{"+ frame dragging",{}}}};
    for(auto& c:cases){c.options.earth_tides=false;c.options.lunar_tides=false;c.options.external_figure_interactions=false;
        c.options.solar_oblateness=false;c.options.solar_radiation_pressure=false;c.options.solar_frame_dragging=false;}
    cases[0].options.earth_degree=2;cases[0].options.moon_degree=2;
    cases[2].options.earth_tides=true;
    cases[3].options.earth_tides=true;cases[3].options.lunar_tides=true;
    for(std::size_t i=4;i<cases.size();++i) {
        cases[i].options.earth_tides=true;cases[i].options.lunar_tides=true;cases[i].options.external_figure_interactions=true;
        if(i>=5)cases[i].options.solar_oblateness=true;
        if(i>=6)cases[i].options.solar_radiation_pressure=true;
        if(i>=7)cases[i].options.solar_frame_dragging=true;
    }
    std::puts("Force ablation, identical seed and +/-90-day window:");
    for(const auto& c:cases) {
        const auto begin=Clock::now();const Ephemeris orbit(kernel,seed,shift(seed,-90),shift(seed,90),{},c.options);
        const double ms=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
        double error=0;
        for(int i=0;i<=720;++i) {
            const auto t=shift(seed,-90+i/4.0);
            error=std::max(error,norm(orbit.icrf(t).position_km-kernel.state(301,399,t).position_km));
        }
        std::printf("  %-20s max %.6f m, %.3f ms build, %zu force calls\n",c.name,error*1000,ms,orbit.stats().force_evaluations);
    }
    // Bound the local predictor's truncation with independent kernel
    // derivatives and true lagged states; no reference states enter Forces.
    double lag_error=0,tidal_error=0;
    for(int i=0;i<16;++i) {
        const auto t=shift(seed,-90+12*i);const auto s=kernel.state(301,399,t);
        constexpr double h=0.01;
        const auto mm=kernel.state(301,399,shift(t,-2*h)).velocity_km_per_day;
        const auto m=kernel.state(301,399,shift(t,-h)).velocity_km_per_day;
        const auto p=kernel.state(301,399,shift(t,h)).velocity_km_per_day;
        const auto pp=kernel.state(301,399,shift(t,2*h)).velocity_km_per_day;
        const auto a=(1/(12*h))*((-1.0)*pp+8*p-8*m+mm);
        const auto jerk=(1/(12*h*h))*((-1.0)*pp+16*p-30*s.velocity_km_per_day+16*m-mm);
        for(double lag:{-0.25,-0.20982243826832672,0.077,0.14519388258636395,0.25}) {
            const auto predicted=lagged_position(s.position_km,s.velocity_km_per_day,a,jerk,lag);
            const auto exact=kernel.state(301,399,shift(t,-lag)).position_km;
            lag_error=std::max(lag_error,norm(predicted-exact));
            for(unsigned order=0;order<3;++order)
                tidal_error=std::max(tidal_error,norm(tide_acceleration(s.position_km,predicted,{0,0,1},order,6378.1366,Central::mu,0.335)
                    -tide_acceleration(s.position_km,exact,{0,0,1},order,6378.1366,Central::mu,0.335)));
        }
    }
    const double projection=0.5*90*90*1000*tidal_error/(398600.43550702266/4902.8001184575496);
    std::printf("Short-lag cubic predictor vs kernel: %.6f m max; tidal-force sensitivity %.9g km/day2 (Earth+Moon GM stress scale); 90-day displacement scale %.9g m\n",lag_error*1000,tidal_error,projection);
    // The predictor is an internal tide-source approximation. Budget its
    // effect on acceleration/position, rather than demand metre accuracy of
    // a source direction hundreds of thousands of kilometres away. 10 um
    // is below this workload's ~50 um integration-refinement residual.
    require(lag_error<0.5 && projection<1e-5,"short-lag predictor outside sampled error budget");
}
void numerical_phases(const char* path) {
    namespace ref=spp::check::phase_reference;
    const auto file=std::filesystem::path(SPP_GRAVITY_REFERENCE).parent_path()/"moon_horizons_phase_tt.csv";
    const auto arcs=ref::read(file.string().c_str());
    const spp::model::spk::Kernel kernel(path);spp::check::ErrorStats errors;
    std::size_t count=0;double build_ms=0,search_ms=0;
    for(int year:{2000,2024,2026,2050}) {
        const auto anchor=spp::Time::from_utc(year,7,1,12);const auto seed=anchor.tdb();
        const auto began=Clock::now();
        const spp::model::lunar_numerical::Ephemeris orbit(kernel,seed,shift(seed,-190),shift(seed,190));
        build_ms+=std::chrono::duration<double,std::milli>(Clock::now()-began).count();
        const auto searching=Clock::now();
        for(const auto& arc:arcs)if(arc.start.to_civil().year==year) {
            const auto actual=ref::event(arc,orbit);
            errors.add(ref::seconds(actual.tt(),arc.root),actual.iso_utc());++count;
        }
        search_ms+=std::chrono::duration<double,std::milli>(Clock::now()-searching).count();
    }
    require(count==62,"numerical phase coverage changed");
    require(errors.report("integrated phase vs DE441 TT","s",0.5),"integrated independent phase gate failed");
    std::printf("Integrated phases: %.3f ms trajectory construction (4 x +/-190 days), %.3f ms searches / %zu events\n",build_ms,search_ms,count);
}
}
int main(int argc,char** argv) {
    try {
        require(argc<=2,"usage: spp_integrator_bench [de440s.bsp]");gravity_controls();numerical();
        if(argc==2){physical(argv[1]);ablation(argv[1]);numerical_phases(argv[1]);}
        std::puts("PASS: numerical propagation accuracy/performance gates");return 0;
    }catch(const std::exception& e){std::fprintf(stderr,"FAIL integration benchmark: %s\n",e.what());return 1;}
}
