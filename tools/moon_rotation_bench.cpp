// Orientation-model, rotating-state and physical-reference gates with timings.
#include "spp/spp.hpp"
#include "spp/model/iau_moon.hpp"
#include "core/vector.hpp"
#include "core/frames.hpp"
#include "core/aberration.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "error_stats.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <format>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
using spp::model::iau_moon::Rotation;
using Epoch=spp::Instant<spp::TDB>;
using namespace spp::detail;
using Clock=std::chrono::steady_clock;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
std::vector<std::vector<double>> read(const std::filesystem::path& path,std::size_t width) {
    std::ifstream file(path);require(file.good(),"cannot open orientation reference");
    std::vector<std::vector<double>> rows;
    for(std::string line;std::getline(file,line);) {
        if(line.empty() || line.front()=='#')continue;
        std::istringstream input(line);std::vector<double> row;
        for(std::string field;std::getline(input,field,',');)row.push_back(std::stod(field));
        require(row.size()==width,"orientation reference column count differs");
        rows.push_back(std::move(row));
    }
    return rows;
}
Epoch epoch(double jd){const auto day=static_cast<std::int64_t>(std::floor(jd));return {day,jd-day};}
Vector direction(double ra,double dec) {
    ra*=spp::deg2rad;dec*=spp::deg2rad;
    return {std::cos(dec)*std::cos(ra),std::cos(dec)*std::sin(ra),std::sin(dec)};
}
void camera_controls() {
    double frame_error=0,tangent_error=0,scale_error=0;
    for(int i=0;i<8;++i) {
        const auto tt=spp::Time::from_utc(-1000+i*575,3,15).tt();
        const auto n=spp::earth::nutation(tt);
        const auto frame=equatorial_frame(tt,n);
        for(Vector v:{Vector{1,0,0},Vector{0,1,0},Vector{0,0,1},Vector{0.3,-0.7,0.2}}) {
            frame_error=std::max(frame_error,norm(frame.to_equatorial(v)-ecliptic_j2000_to_equatorial(v,tt,n)));
            frame_error=std::max(frame_error,norm(frame.to_ecliptic(v)-equatorial_to_ecliptic_j2000(v,tt,n)));
        }
    }
    for(int i=0;i<32;++i)for(double beta:{0.0,1e-4,0.01,0.6}) {
        const double l=i*spp::two_pi/32;
        const Vector n{std::cos(l),std::sin(l),0},tangent{-std::sin(l),std::cos(l),0};
        const Vector v=light_speed_km_day*Vector{0.6*beta,0,-0.8*beta};
        const auto apparent=aberration(n,v),mapped=aberration_tangent(n,apparent,tangent,v);
        constexpr double h=1e-6;
        const auto difference=(0.5/h)*(aberration(n+h*tangent,v)-aberration(n-h*tangent,v));
        tangent_error=std::max(tangent_error,norm(mapped-difference));
        const double scale=std::sqrt(1-beta*beta)/(1+beta*0.6*n[0]);
        scale_error=std::max(scale_error,std::abs(norm(mapped)-scale));
        require(std::abs(dot(mapped,apparent))<1e-14,"camera tangent not in apparent sky plane");
    }
    std::printf("Observer basis controls: frame reuse %.9g, aberration tangent %.9g, conformal scale %.9g\n",frame_error,tangent_error,scale_error);
    require(frame_error<4e-15 && tangent_error<1e-8 && scale_error<1e-14,"observer basis numerical gate failed");
}
spp::moon::BodyRotation recover_rotation(const std::vector<double>& v) {
    const Vector z=direction(v[6],v[7]),observer=(-1.0)*direction(v[2],v[3]);
    const auto projected=observer-dot(observer,z)*z;
    const Vector u=(1.0/norm(projected))*projected,w=cross(z,u);
    const double lon=v[4]*spp::deg2rad;
    spp::moon::BodyRotation result;
    result.matrix={std::cos(lon)*u-std::sin(lon)*w,std::sin(lon)*u+std::cos(lon)*w,z};
    return result;
}
void frame_reference(const std::filesystem::path& path) {
    const auto bias=read(path,3);require(bias.size()==3,"ERFA bias reference changed");
    const double eps=84381.448/3600*spp::deg2rad,c=std::cos(eps),s=std::sin(eps);
    const spp::RotationMatrix tilt{{{1,0,0},{0,c,-s},{0,s,c}}};
    const auto actual=spp::frames::mean_ecliptic_j2000_to_icrf();
    double error=0;
    for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
        double expected=0;
        for(std::size_t k=0;k<3;++k)expected+=bias[k][i]*tilt[k][j];
        error=std::max(error,std::abs(actual.matrix[i][j]-expected));
        require(actual.derivative_per_day[i][j]==0,"fixed J2000 adapter has a spurious rate");
    }
    std::printf("J2000 ecliptic -> ICRF adapter vs independent ERFA frame bias: %.9g matrix error\n",error);
    require(error<2e-15,"ICRF adapter has the wrong frame-bias direction");
}
void compatibility(const char* path) {
    const Rotation model;
    static_assert(!std::is_invocable_v<Rotation,spp::Instant<spp::TT>>);
    using R=spp::FrameRotation<spp::frame::Icrf,spp::frame::MoonMeanEarth>;
    static_assert(!std::is_invocable_v<R,spp::GeocentricEclipticJ2000State>);
    double max_angle=0,max_matrix=0,max_rate=0,max_roundtrip=0,max_velocity=0;
    const auto rows=read(path,13);require(rows.size()==8,"IAU reference coverage changed");
    for(const auto& v:rows) {
        const auto t=epoch(v[0]);const auto o=model(t);
        max_angle=std::max({max_angle,std::abs(spp::wrap180(o.pole_ra_deg-v[1])),
            std::abs(o.pole_dec_deg-v[2]),std::abs(spp::wrap180(o.prime_meridian_deg-v[3]))});
        const auto& r=o.rotation;
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
            max_matrix=std::max(max_matrix,std::abs(r.matrix[i][j]-v[4+3*i+j]));
        const spp::CartesianState<spp::origin::Geocentric,spp::frame::Icrf> state{{350000,-200000,10000},{40000,70000,-1000}};
        const auto restored=r.inverse()(r(state));
        max_roundtrip=std::max(max_roundtrip,norm(restored.position_km-state.position_km));
        require(norm(restored.velocity_km_per_day-state.velocity_km_per_day)<1e-8,"rotating-frame state inverse lost velocity terms");
        const auto bridge=spp::frames::mean_ecliptic_j2000_to_icrf();
        const spp::GeocentricEclipticJ2000State ecliptic{state.position_km,state.velocity_km_per_day};
        const auto composed=(r*bridge)(ecliptic),sequential=r(bridge(ecliptic));
        require(norm(composed.position_km-sequential.position_km)<1e-9
            && norm(composed.velocity_km_per_day-sequential.velocity_km_per_day)<1e-8,
            "typed ecliptic/ICRF/body rotation composition failed");
        const auto identity=r.inverse()*r;
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
            require(std::abs(identity.matrix[i][j]-(i==j?1:0))<1e-14,"rotation composition not identity");
            require(std::abs(identity.derivative_per_day[i][j])<1e-14,"rotation derivative composition not zero");
        }
        if(t.day>Rotation::first_epoch.day && t.day<Rotation::last_epoch.day) {
            constexpr double h=0.001;
            const auto mm=model({t.day,t.frac-2*h}).rotation,m=model({t.day,t.frac-h}).rotation;
            const auto p=model({t.day,t.frac+h}).rotation,pp=model({t.day,t.frac+2*h}).rotation;
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
                const double rate=(-pp.matrix[i][j]+8*p.matrix[i][j]-8*m.matrix[i][j]+mm.matrix[i][j])/(12*h);
                max_rate=std::max(max_rate,std::abs(rate-r.derivative_per_day[i][j]));
            }
            auto moving_point=[&](double offset) {
                const spp::CartesianPosition<spp::origin::Geocentric,spp::frame::Icrf> point{
                    state.position_km+offset*state.velocity_km_per_day};
                return model({t.day,t.frac+offset}).rotation(point).position_km;
            };
            const Vector velocity=(1.0/(12*h))*((-1.0)*moving_point(2*h)+8*moving_point(h)
                -8*moving_point(-h)+moving_point(-2*h));
            max_velocity=std::max(max_velocity,norm(velocity-r(state).velocity_km_per_day));
        }
    }
    std::printf("IAU/PCK compatibility: max angle %.9g deg, matrix %.9g, rate %.9g/day, state roundtrip %.9g m, velocity derivative %.9g km/day\n",
        max_angle,max_matrix,max_rate,max_roundtrip*1000,max_velocity);
    require(max_angle<1e-8 && max_matrix<5e-10 && max_rate<1e-8 && max_roundtrip<1e-9
        && max_velocity<1e-5,"IAU rotation numerical gate failed");
}
void physical(const char* path) {
    const Rotation model;spp::check::ErrorStats axes;
    double latitude_control=0,pole_error=0;
    const auto rows=read(path,8);require(rows.size()==7,"JPL rotation coverage changed");
    for(const auto& v:rows) {
        const auto day=static_cast<std::int64_t>(std::floor(v[0]));
        const auto t=spp::tdb_from_tt({day,v[0]-day-v[1]/1440});
        const auto r=model(t).rotation;
        const Vector z=direction(v[6],v[7]),observer=(-1.0)*direction(v[2],v[3]);
        latitude_control=std::max(latitude_control,std::abs(std::asin(std::clamp(dot(z,observer),-1.0,1.0))*spp::rad2deg-v[5]));
        double squared=0;
        const auto reference=recover_rotation(v).matrix;
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
            const double d=r.matrix[i][j]-reference[i][j];squared+=d*d;
        }
        pole_error=std::max(pole_error,angle(r.matrix[2],z)*spp::rad2deg);
        const double error=2*std::asin(std::min(1.0,std::sqrt(squared/8)))*spp::rad2deg;
        axes.add(error,std::format("TT reception JD {:.8f}",v[0]));
    }
    std::printf("JPL orientation reference control: sub-observer latitude %.9g arcsec; max IAU pole error %.6f arcsec\n",
        latitude_control*3600,pole_error*3600);
    require(latitude_control<0.001,"JPL pole/sub-observer/astrometric frame control failed");
    require(axes.report("IAU vs DE441 ME axes","deg",0.1),"IAU rotation outside sampled physical budget");
}
void surface_cases(const std::filesystem::path& surface_path,const std::filesystem::path& rotation_path) {
    struct Case {spp::Time t;std::optional<spp::Observer> observer;spp::moon::Surface expected;};
    struct Frozen {
        spp::moon::BodyRotation matrix;int calls=0;Epoch last{};
        const spp::moon::BodyRotation& operator()(Epoch t){++calls;last=t;return matrix;}
    };
    const auto references=read(rotation_path,8);
    spp::model::elp_mpp02::CachedEphemeris model;
    std::ifstream file(surface_path);require(file.good(),"cannot open surface reference");
    std::vector<Case> cases;spp::check::ErrorStats observer_error,solar_error,pole_error,geometry;
    for(std::string line;std::getline(file,line);) {
        if(line.empty() || line.front()=='#')continue;
        std::istringstream input(line);std::string source;std::getline(input,source,',');
        std::vector<double> v;
        for(std::string value;std::getline(input,value,',');)v.push_back(std::stod(value));
        require(v.size()==10,"surface reference column count differs");
        const auto day=static_cast<std::int64_t>(std::floor(v[0]));
        const bool geo=source=="geocenter";
        auto t=spp::Time::from_jd(day,v[0]-day,geo?std::nullopt:std::optional<double>{v[4]});
        for(int i=0;i<4;++i) {const auto tt=t.tt();t=t.shifted(static_cast<double>(day-tt.day)+(v[0]-day-tt.frac));}
        const std::optional<spp::Observer> obs=geo?std::nullopt:std::optional<spp::Observer>{{v[1],v[2],v[3]}};
        auto actual=geo?spp::moon::surface(t,model):spp::moon::surface(t,*obs,model);
        const auto label=std::format("{} TT JD {:.8f}",source,v[0]);
        observer_error.add(std::hypot(spp::wrap180(actual.subobserver_lon_deg-v[5]),actual.subobserver_lat_deg-v[6]),label);
        solar_error.add(std::hypot(spp::wrap180(actual.subsolar_lon_deg-v[7]),actual.subsolar_lat_deg-v[8]),label);
        require(actual.pole_defined,"Moon pole unexpectedly undefined");
        pole_error.add(spp::wrap180(actual.pole_position_angle_deg-v[9]),label);
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
            require(std::abs(dot(actual.body_to_view[i],actual.body_to_view[j])-(i==j?1:0))<1e-12,"body view matrix not orthonormal");
        const auto observed=geo?spp::moon::Surface{}:spp::moon::observe(t,*obs,model).surface;
        if(!geo)require(observed.body_to_view==actual.body_to_view
            && observed.subobserver_lon_deg==actual.subobserver_lon_deg,"standalone/complete surface paths differ");
        const auto reference=std::find_if(references.begin(),references.end(),[&](const auto& row){return row[0]==v[0];});
        require(reference!=references.end(),"surface oracle has no reference rotation");
        Frozen spin{recover_rotation(*reference)};
        const auto fixed=geo?spp::moon::surface(t,model,spin):spp::moon::surface(t,*obs,model,spin);
        require(spin.calls==1,"rotation provider copied or evaluated repeatedly");
        const auto reception=t.tdb();
        const double delay=(static_cast<double>(reception.day-spin.last.day)+reception.frac-spin.last.frac)*spp::seconds_per_day;
        require(delay>1 && delay<1.5,"rotation provider not evaluated at lunar reflection time");
        const auto body_observer=direction(fixed.subobserver_lon_deg,fixed.subobserver_lat_deg);
        require(norm(fixed.body_to_view[2]-body_observer)<1e-13,"surface camera moved the sub-observer point");
        if(!geo) {
            const auto combined=spp::moon::observe(t,*obs,model,spin);
            require(spin.calls==2 && combined.surface.body_to_view==fixed.body_to_view,
                "complete observation ignored the explicit rotation provider");
            for(int invalid=0;invalid<3;++invalid) {
                auto bad=spin;
                if(invalid==0)bad.matrix.matrix[0][0]+=0.1;
                if(invalid==1)bad.matrix.matrix[0]=(-1.0)*bad.matrix.matrix[0];
                if(invalid==2)bad.matrix.derivative_per_day[0][0]=std::numeric_limits<double>::infinity();
                bool rejected=false;
                try{(void)spp::moon::surface(t,*obs,model,bad);}catch(const std::invalid_argument&){rejected=true;}
                require(rejected,"invalid rotation provider accepted");
            }
        }
        if(geo) {
            geometry.add(3600*std::hypot(spp::wrap180(fixed.subobserver_lon_deg-v[5]),fixed.subobserver_lat_deg-v[6]),label+" observer");
            geometry.add(3600*std::hypot(spp::wrap180(fixed.subsolar_lon_deg-v[7]),fixed.subsolar_lat_deg-v[8]),label+" Sun");
        }
        cases.push_back({t,obs,actual});
    }
    require(cases.size()==13,"surface coverage changed");
    const bool observer_ok=observer_error.report("sub-observer vs DE441","deg",0.01);
    const bool solar_ok=solar_error.report("subsolar vs DE441","deg",0.01);
    const bool pole_ok=pole_error.report("pole bearing vs DE441","deg",0.02);
    const bool geometry_ok=geometry.report("frozen reference geometry","arcsec",0.5);
    require(observer_ok && solar_ok && pole_ok && geometry_ok,"surface/reference angular gate failed");
    const auto calls=model.source_evaluations();double checksum=0;
    const auto start=Clock::now();
    for(int repeat=0;repeat<128;++repeat)for(const auto& q:cases) {
        const auto s=q.observer?spp::moon::surface(q.t,*q.observer,model):spp::moon::surface(q.t,model);
        require(s.body_to_view==q.expected.body_to_view && s.subobserver_lon_deg==q.expected.subobserver_lon_deg,
            "surface repeat changed its result");
        checksum+=s.subobserver_lon_deg+s.subobserver_lat_deg+s.subsolar_lon_deg+s.subsolar_lat_deg+s.pole_position_angle_deg;
        for(const auto& row:s.body_to_view)for(double x:row)checksum+=x;
    }
    const double us=std::chrono::duration<double,std::micro>(Clock::now()-start).count()/(128*cases.size());
    std::printf("Surface + view orientation: %.3f us/query resident; %llu ELP refills; checksum %.9f\n",
        us,static_cast<unsigned long long>(model.source_evaluations()-calls),checksum);
    require(std::isfinite(checksum) && model.source_evaluations()==calls,"surface cache workload unnecessarily refilled");
}
void speed() {
    const Rotation model;std::vector<Epoch> input;
    for(int i=0;i<4096;++i)input.push_back({2451545+(i*1543)%36525,(i%73)/73.0});
    double checksum=0;std::vector<double> times;
    for(int batch=0;batch<9;++batch) {
        const auto begin=Clock::now();
        for(const auto t:input) {
            const auto o=model(t);
            checksum+=o.pole_ra_deg+o.pole_dec_deg+o.prime_meridian_deg
                +o.pole_ra_rate_deg_day+o.pole_dec_rate_deg_day+o.prime_meridian_rate_deg_day;
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
                checksum+=o.rotation.matrix[i][j]+o.rotation.derivative_per_day[i][j];
        }
        times.push_back(std::chrono::duration<double,std::micro>(Clock::now()-begin).count()/input.size());
    }
    std::sort(times.begin(),times.end());require(std::isfinite(checksum),"invalid rotation benchmark checksum");
    std::printf("IAU rotation + derivative: median %.3f us/query, min %.3f, max %.3f; 4096 dates x 9 batches; checksum %.9f\n",
        times[4],times.front(),times.back(),checksum);
}
}
int main(int argc,char** argv) {
    try {
        require(argc==3,"usage: spp_moon_rotation_bench moon_iau_rotation.csv moon_rotation_horizons.csv");
        frame_reference(std::filesystem::path(argv[1]).parent_path()/"frame_bias_erfa.csv");
        camera_controls();compatibility(argv[1]);physical(argv[2]);
        surface_cases(std::filesystem::path(argv[1]).parent_path()/"moon_surface_horizons.csv",argv[2]);speed();
        std::puts("PASS: lunar rotation model, frame-rate and performance gates");return 0;
    } catch(const std::exception& e){std::fprintf(stderr,"FAIL rotation benchmark: %s\n",e.what());return 1;}
}
