#include "spp/model/lunar_numerical.hpp"
#include "spp/earth.hpp"
#include "spp/model/iau_moon.hpp"
#include "../core/vector.hpp"
#include "../core/frames.hpp"
#include "lunar_gravity.hpp"
#include "../core/tides.hpp"
#include "../core/chebyshev.hpp"
#include "spp/core/detail/split.hpp"
#include "solar_effects.hpp"
#include <cmath>
#include <stdexcept>

namespace spp::model::lunar_numerical {
namespace {
using namespace spp::detail;
constexpr double day2=seconds_per_day*seconds_per_day;
// NAIF gm_de440.tpc, km^3/s^2. Earth, Moon, Sun, and system barycenters.
constexpr std::array<int,11> ids{399,301,10,1,2,4,5,6,7,8,9};
constexpr std::array<double,11> gm{398600.43550702266,4902.8001184575496,132712440041.27942,
    22031.868551400003,324858.592,42828.375815756102,126712764.09999998,
    37940584.841799997,5794556.3999999985,6836527.1005803989,975.5};
using States=std::array<State,11>;
using Vectors=std::array<Vector,11>;
// Park et al. (2021), equation 27, beta=gamma=1. This returns only 1PN
// corrections: Newtonian gravity is accumulated separately and stably.
Vector pn(std::size_t i,const States& states,const Vectors& accelerations,const std::array<double,11>& potential,
    const std::array<double,11>& inverse) {
    constexpr double invc2=1/(299792.458*seconds_per_day*299792.458*seconds_per_day);
    const auto& ri=states[i].position_km;const auto& vi=states[i].velocity_km_per_day;
    Vector result{};
    for(std::size_t j=0;j<states.size();++j)if(i!=j) {
        const auto& vj=states[j].velocity_km_per_day;const auto& aj=accelerations[j];
        const Vector dr=states[j].position_km-ri;const double inv=inverse[j],mu=gm[j]*day2;
        const double projection=dot(dr,vj)*inv;
        const double bracket=-4*potential[i]-potential[j]+dot(vi,vi)+2*dot(vj,vj)-4*dot(vi,vj)
            -1.5*projection*projection+0.5*dot(dr,aj);
        const double scale=mu*inv*inv*inv*invc2;
        result=result+(scale*bracket)*dr;
        result=result+(-scale*dot(dr,4*vi-3*vj))*(vi-vj);
        result=result+(3.5*mu*inv*invc2)*aj;
    }
    return result;
}
State seed_state(const spk::Kernel& kernel,Instant<TDB> seed) {
    const auto s=kernel.state(301,399,seed);return {s.position_km,s.velocity_km_per_day};
}
Instant<TT> terrestrial_time(Instant<TDB> t) {
    // Invert Spp's periodic TDB-TT map; using the same numerical fields as
    // an initial guess is not a relabeling of the final physical instant.
    Instant<TT> tt{t.day,t.frac};
    for(int i=0;i<2;++i) {
        const auto estimate=tdb_from_tt(tt);
        tt.frac+=static_cast<double>(t.day-estimate.day)+(t.frac-estimate.frac);
    }
    const double whole=std::floor(tt.frac);return {tt.day+static_cast<std::int64_t>(whole),tt.frac-whole};
}
struct PrincipalAxes {};
using LunarRotation=FrameRotation<frame::Icrf,PrincipalAxes>;
LunarRotation lunar_axes(Instant<TDB> t) {
    const auto me=iau_moon::Rotation{}(t).rotation;
    // NAIF moon_de440_250416.tf: ME_DE440_ME421 -> PA_DE440.
    // Passive R3(67.8526") R2(78.6944") R1(0.2785").
    const double a=67.8526*deg2rad/3600,b=78.6944*deg2rad/3600,c=0.2785*deg2rad/3600;
    LunarRotation result{};
    for(std::size_t j=0;j<3;++j) {
        for(bool rate:{false,true}) {
            const auto& input=rate?me.derivative_per_day:me.matrix;
            auto& output=rate?result.derivative_per_day:result.matrix;
            auto v=rotate_x({input[0][j],input[1][j],input[2][j]},-c);
            v={std::cos(b)*v[0]-std::sin(b)*v[2],v[1],std::sin(b)*v[0]+std::cos(b)*v[2]};
            v=rotate_z(v,-a);
            for(std::size_t i=0;i<3;++i)output[i][j]=v[i];
        }
    }
    return result;
}
}
// Only prescribed, time-dependent inputs live here. No propagated state,
// acceleration, lunar kernel state, or state-dependent PN/tidal term is fitted.
struct Forces::Environment {
    struct Value {
        std::array<spk::State,10> external{}; // EMB/SSB followed by planets/EMB
        std::array<Vector,3> delayed_sun{};
        Vector pole{};
        RotationMatrix axes{};
        LunarRotation delayed_axes{};
        std::array<Vector,9> external_acceleration{};
        std::array<double,9> external_potential{};
    };
    using Series=ChebyshevSeries<135,12>;
    struct Segment {Instant<TDB> first;double width;bool cacheable;Series series;};
    Instant<TDB> first,last;
    std::vector<Segment> segments;
    static constexpr double interval_days=4;
    static double difference(Instant<TDB> a,Instant<TDB> b){return static_cast<double>(a.day-b.day)+(a.frac-b.frac);}
    static Value sample(const Forces& forces,Instant<TDB> t,bool summarize=false) {
        Value value;const auto& options=forces.options_;
        constexpr std::array<int,9> targets{10,1,2,4,5,6,7,8,9};
        forces.kernel_.states(targets,3,t,std::span<spk::State>(value.external).subspan(1));
        if(options.relativity)value.external[0]=forces.kernel_.state(3,0,t);
        if(options.earth_oblateness || options.earth_tides) {
            const auto tt=terrestrial_time(t);const auto nut=earth::nutation(tt);
            const auto pole0=equatorial_frame(tt,nut).to_ecliptic({0,0,1});
            using Direction=CartesianPosition<origin::Geocentric,frame::MeanEclipticJ2000>;
            value.pole=frames::mean_ecliptic_j2000_to_icrf()(Direction{pole0}).position_km;
        }
        if(options.earth_tides)for(unsigned j=0;j<3;++j)
            value.delayed_sun[j]=forces.kernel_.state(10,3,{t.day,t.frac-options.earth_tide.orbit_lag_days[j]}).position_km;
        if(options.lunar_figure) {
            value.axes=lunar_axes(t).matrix;
            if(options.lunar_tides)value.delayed_axes=lunar_axes({t.day,t.frac-options.moon_lag_days});
        }
        if(options.relativity && summarize)for(std::size_t i=0;i<9;++i)for(std::size_t j=i+1;j<9;++j) {
            const auto d=value.external[j+1].position_km-value.external[i+1].position_km;
            const double inv=1/norm(d),inv3=inv*inv*inv;
            value.external_acceleration[i]=value.external_acceleration[i]+(gm[j+2]*day2*inv3)*d;
            value.external_acceleration[j]=value.external_acceleration[j]-(gm[i+2]*day2*inv3)*d;
            value.external_potential[i]+=gm[j+2]*day2*inv;value.external_potential[j]+=gm[i+2]*day2*inv;
        }
        return value;
    }
    static Series::Value pack(const Value& v) {
        Series::Value out;std::size_t n=0;
        auto put=[&](const Vector& a){for(double x:a)out[n++]=x;};
        for(const auto& s:v.external){put(s.position_km);put(s.velocity_km_per_day);}
        for(const auto& s:v.delayed_sun)put(s);
        put(v.pole);for(const auto& s:v.axes)put(s);
        for(const auto& s:v.delayed_axes.matrix)put(s);
        for(const auto& s:v.delayed_axes.derivative_per_day)put(s);
        for(const auto& s:v.external_acceleration)put(s);
        for(double x:v.external_potential)out[n++]=x;
        return out;
    }
    static Value unpack(const Series::Value& in) {
        Value out;std::size_t n=0;
        auto get=[&](Vector& a){for(double& x:a)x=in[n++];};
        for(auto& s:out.external){get(s.position_km);get(s.velocity_km_per_day);}
        for(auto& s:out.delayed_sun)get(s);
        get(out.pole);for(auto& s:out.axes)get(s);
        for(auto& s:out.delayed_axes.matrix)get(s);
        for(auto& s:out.delayed_axes.derivative_per_day)get(s);
        for(auto& s:out.external_acceleration)get(s);
        for(double& x:out.external_potential)x=in[n++];
        return out;
    }
    Environment(const Forces& forces,Instant<TDB> begin,Instant<TDB> end):first(begin),last(end) {
        // Chebyshev nodes are interior. Explicit endpoint samples preserve
        // kernel/rotation coverage errors at the actual requested boundaries.
        (void)sample(forces,first);(void)sample(forces,last);
        const double span=difference(last,first);
        const auto count=static_cast<std::size_t>(std::ceil(span/interval_days));segments.reserve(count);
        for(std::size_t i=0;i<count;++i) {
            const double offset=i*interval_days,width=std::min(interval_days,span-offset);
            const auto start=spp::detail::normalized(Instant<TDB>{first.day,first.frac+offset});
            bool cacheable=true;
            // A new segment can change SPK priority discontinuously. Do not
            // smooth across any such boundary, including delayed solar input.
            std::array<double,4> lags{};
            if(forces.options_.earth_tides)for(unsigned j=0;j<3;++j)lags[j+1]=forces.options_.earth_tide.orbit_lag_days[j];
            for(const auto& segment:forces.kernel_.segments())for(const auto boundary:{segment.first,segment.last})for(double lag:lags) {
                const double at=difference(boundary,start)+lag;
                if(at>=0 && at<=width)cacheable=false;
            }
            Series series;
            if(cacheable)series=Series::fit([&](double x){
                return pack(sample(forces,{start.day,start.frac+0.5*width*(x+1)},true));
            });
            segments.push_back({start,width,cacheable,series});
        }
    }
    const Segment& segment(Instant<TDB> t)const {
        const double offset=difference(t,first);
        const auto index=std::min(static_cast<std::size_t>(offset/interval_days),segments.size()-1);
        return segments[index];
    }
    Value at(Instant<TDB> t)const {
        const auto& s=segment(t);
        return unpack(s.series(2*difference(t,s.first)/s.width-1));
    }
};
Forces Forces::prepared(Instant<TDB> first,Instant<TDB> last)const {
    Forces result=*this;result.environment_.reset();
    if(!options_.environment_cache)return result;
    for(const auto t:{first,last})if(!std::isfinite(t.frac) || std::abs(t.frac)>1e9 || t.day<-1000000000 || t.day>1000000000)
        throw std::invalid_argument("spp lunar forces: unusable environment domain");
    first=spp::detail::normalized(first);last=spp::detail::normalized(last);
    const double span=Environment::difference(last,first);
    if(span<0)throw std::invalid_argument("spp lunar forces: reversed environment domain");
    // Bound auxiliary storage to ~62 MB. Very long research integrations use
    // the exact source instead; no extrapolation or unbounded allocation.
    if(span>0 && span<=4096*Environment::interval_days)
        result.environment_=std::make_shared<Environment>(*this,first,last);
    return result;
}
std::size_t Forces::environment_bytes()const {
    return environment_?sizeof(Environment)+environment_->segments.capacity()*sizeof(Environment::Segment):0;
}
std::size_t Forces::environment_samples()const {
    std::size_t count=0;
    if(environment_) {
        count=2;
        for(const auto& s:environment_->segments)if(s.cacheable)count+=Environment::Series::nodes;
    }
    return count;
}
Forces::Forces(spk::Kernel kernel,ForceOptions options):kernel_(std::move(kernel)),options_(options) {
    if(!std::isfinite(options_.earth_j2) || std::abs(options_.earth_j2)>0.1
        || !std::isfinite(options_.earth_radius_km) || options_.earth_radius_km<=0
        || !std::isfinite(options_.moon_radius_km) || options_.moon_radius_km<=0
        || !std::isfinite(options_.moon_j2) || std::abs(options_.moon_j2)>0.1
        || !std::isfinite(options_.moon_c22) || std::abs(options_.moon_c22)>0.1)
        throw std::invalid_argument("spp lunar forces: invalid body figure");
    if(options_.earth_degree<2 || options_.earth_degree>5 || options_.moon_degree<2 || options_.moon_degree>6)
        throw std::invalid_argument("spp lunar forces: invalid gravity degree");
    for(unsigned j=0;j<3;++j) {
        const auto& tide=options_.earth_tide;
        if(!std::isfinite(tide.love[j]) || tide.love[j]<0 || tide.love[j]>1
            || !std::isfinite(tide.orbit_lag_days[j]) || std::abs(tide.orbit_lag_days[j])>0.25
            || !std::isfinite(tide.rotation_lag_days[j]) || std::abs(tide.rotation_lag_days[j])>0.1)
            throw std::invalid_argument("spp lunar forces: invalid Earth tidal response or short-lag domain");
    }
    if(options_.earth_tide.rotation_lag_days[0]!=0)throw std::invalid_argument("spp lunar forces: zonal tide has no rotational lag");
    if(!std::isfinite(options_.moon_love) || options_.moon_love<0 || options_.moon_love>1
        || !std::isfinite(options_.moon_lag_days) || std::abs(options_.moon_lag_days)>0.25
        || !std::isfinite(options_.moon_reference_spin_rad_day) || options_.moon_reference_spin_rad_day<=0)
        throw std::invalid_argument("spp lunar forces: invalid lunar tidal response");
    for(double fraction:{options_.earth_srp_fraction,options_.moon_srp_fraction})
        if(!std::isfinite(fraction) || fraction<0 || fraction>1)
            throw std::invalid_argument("spp lunar forces: invalid radiation pressure fraction");
}
Vector Forces::operator()(Instant<TDB> t,const State& state)const {
    if(!std::isfinite(t.frac))throw std::invalid_argument("spp lunar forces: non-finite TDB date");
    if(t.day<-1000000000 || t.day>1000000000 || std::abs(t.frac)>1e9)
        throw std::out_of_range("spp lunar forces: unusable TDB date");
    const Vector r=state.position_km;const double distance=norm(r);
    if(!std::isfinite(distance) || distance<options_.earth_radius_km)
        throw std::invalid_argument("spp lunar forces: invalid/colliding Moon position");
    for(double v:state.velocity_km_per_day)if(!std::isfinite(v))throw std::invalid_argument("spp lunar forces: invalid velocity");
    Vector result=(-(gm[0]+gm[1])*day2/(distance*distance*distance))*r;
    const Vector v=state.velocity_km_per_day;
    Vector jerk=(-(gm[0]+gm[1])*day2/(distance*distance*distance))*(v-(3*dot(r,v)/(distance*distance))*r);
    States bodies{};
    Vectors relative_positions{};
    Vector sun_velocity{};
    // Prescribe the Earth-Moon barycenter, then derive both members from
    // the propagated separation. The kernel's lunar recoil in its Earth
    // segment must not feed the reference Moon trajectory back into forces.
    constexpr double fraction=gm[1]/(gm[0]+gm[1]);
    const bool resident=environment_ && Environment::difference(t,environment_->first)>=0
        && Environment::difference(environment_->last,t)>=0 && environment_->segment(t).cacheable;
    const auto environment=resident?environment_->at(t):Environment::sample(*this,t);
    const auto& emb=environment.external[0];
    const State earth{emb.position_km-fraction*r,emb.velocity_km_per_day-fraction*state.velocity_km_per_day};
    bodies[0]=earth;
    bodies[1]={emb.position_km+(1-fraction)*r,emb.velocity_km_per_day+(1-fraction)*state.velocity_km_per_day};
    for(std::size_t i=2;i<ids.size();++i) {
        const auto& relative=environment.external[i-1];
        const Vector q=relative.position_km+fraction*r,d=q-r;
        if((options_.external_figure_interactions && i<=7) || i==2)relative_positions[i]=q;
        if(i==2)sun_velocity=relative.velocity_km_per_day+fraction*v;
        const double q2=dot(q,q),d2=dot(d,d),d3=d2*std::sqrt(d2);
        // Rationalize (1+delta)^(3/2)-1 to preserve the small differential
        // acceleration without a logarithm/exponential pair for every body.
        const double delta=(dot(r,r)-2*dot(q,r))/q2;
        if(delta<=-1 || d2<=0)throw std::invalid_argument("spp lunar forces: perturbing-body collision");
        const double f=delta*(3+delta*(3+delta))/(1+(1+delta)*std::sqrt(1+delta));
        result=result+(-gm[i]*day2/d3)*(r+f*q);
        if(options_.earth_tides || options_.lunar_tides) {
            const Vector qv=relative.velocity_km_per_day+fraction*v,dv=qv-v;
            jerk=jerk+(gm[i]*day2)*((1/d3)*(dv-(3*dot(d,dv)/d2)*d)
                -(1/(q2*std::sqrt(q2)))*(qv-(3*dot(q,qv)/q2)*q));
        }
        if(options_.relativity)bodies[i]={emb.position_km+relative.position_km,emb.velocity_km_per_day+relative.velocity_km_per_day};
    }
    const auto newtonian=result;
    const Vector sun_to_earth=(-1.0)*relative_positions[2],sun_to_moon=sun_to_earth+r;
    if(options_.solar_oblateness)result=result+de440_solar::oblateness(sun_to_moon,gm[2]*day2)
        -de440_solar::oblateness(sun_to_earth,gm[2]*day2);
    if(options_.solar_radiation_pressure)result=result+de440_solar::radiation(sun_to_moon,gm[2]*day2,options_.moon_srp_fraction)
        -de440_solar::radiation(sun_to_earth,gm[2]*day2,options_.earth_srp_fraction);
    if(options_.relativity && options_.solar_frame_dragging)
        result=result+de440_solar::frame_dragging(sun_to_moon,v-sun_velocity,gm[2]*day2)
            -de440_solar::frame_dragging(sun_to_earth,(-1.0)*sun_velocity,gm[2]*day2);
    if(options_.earth_oblateness || options_.earth_tides) {
        // True pole of date, with explicit ICRF/J2000 bias adapter.
        // Zonal gravity is axisymmetric; only tidal lag needs the spin rate.
        const Vector pole=environment.pole;
        if(options_.earth_oblateness) {
            const std::array<double,6> j{0,0,options_.earth_j2,-2.53241e-6,-1.619898e-6,-2.277345e-7};
            result=result+zonal_acceleration(r,pole,options_.earth_radius_km,(gm[0]+gm[1])*day2,options_.earth_degree,j);
            if(options_.external_figure_interactions)for(std::size_t i=2;i<=7;++i)
                result=result+zonal_acceleration(relative_positions[i],pole,options_.earth_radius_km,
                    gm[i]*day2,options_.earth_degree,j);
        }
        if(options_.earth_tides)for(unsigned j=0;j<3;++j) {
            const auto& tide=options_.earth_tide;
            const double lag=tide.orbit_lag_days[j];
            const auto moon=lagged_position(r,v,newtonian,jerk,lag);
            const Vector sun=environment.delayed_sun[j]+fraction*moon;
            // DE430 eqs 13/31: Rz(-omega*tau) is a passive matrix,
            // hence the corresponding vector rotates by +omega*tau.
            const double angle=earth_rotation_rad_s*seconds_per_day*tide.rotation_lag_days[j];
            const double reaction=1+gm[1]/gm[0];
            result=result+tide_acceleration(r,rotate_about(moon,pole,angle),pole,j,
                options_.earth_radius_km,reaction*gm[1]*day2,tide.love[j]);
            result=result+tide_acceleration(r,rotate_about(sun,pole,angle),pole,j,
                options_.earth_radius_km,reaction*gm[2]*day2,tide.love[j]);
        }
    }
    if(options_.relativity) {
        Vectors accelerations{};std::array<double,11> potential{};
        std::array<std::array<double,11>,2> inverse{};
        for(std::size_t i=0;i<(resident?2:bodies.size());++i)for(std::size_t j=i+1;j<bodies.size();++j) {
            const auto d=bodies[j].position_km-bodies[i].position_km;const double inv=1/norm(d),inv3=inv*inv*inv;
            if(i<2){inverse[i][j]=inv;if(j<2)inverse[j][i]=inv;}
            accelerations[i]=accelerations[i]+(gm[j]*day2*inv3)*d;
            accelerations[j]=accelerations[j]-(gm[i]*day2*inv3)*d;
            potential[i]+=gm[j]*day2*inv;potential[j]+=gm[i]*day2*inv;
        }
        if(resident)for(std::size_t i=2;i<bodies.size();++i) {
            accelerations[i]=accelerations[i]+environment.external_acceleration[i-2];
            potential[i]+=environment.external_potential[i-2];
        }
        result=result+pn(1,bodies,accelerations,potential,inverse[1])-pn(0,bodies,accelerations,potential,inverse[0]);
    }
    if(options_.lunar_figure) {
        const auto& axes=environment.axes;
        const Vector body{dot(axes[0],r),dot(axes[1],r),dot(axes[2],r)};
        auto c=detail::moon_cosine;auto s=detail::moon_sine;
        c[2][0]=-options_.moon_j2;c[2][2]=options_.moon_c22;
        if(options_.lunar_tides) {
            const double lag=options_.moon_lag_days;
            const auto& delayed=environment.delayed_axes;
            const auto earth_direction=(-1.0)*lagged_position(r,v,newtonian,jerk,lag);
            const Vector q{dot(delayed.matrix[0],earth_direction),dot(delayed.matrix[1],earth_direction),dot(delayed.matrix[2],earth_direction)};
            const double qnorm=norm(q);
            const Vector omega{dot(delayed.derivative_per_day[1],delayed.matrix[2]),
                dot(delayed.derivative_per_day[2],delayed.matrix[0]),dot(delayed.derivative_per_day[0],delayed.matrix[1])};
            const double r3=options_.moon_radius_km*options_.moon_radius_km*options_.moon_radius_km,
                n2=options_.moon_reference_spin_rad_day*options_.moon_reference_spin_rad_day;
            const double q2=qnorm*qnorm;
            const double tidal=-options_.moon_love*(gm[0]/gm[1])*r3/(q2*q2*qnorm);
            const double spin=options_.moon_love*r3/(3*gm[1]*day2);
            RotationMatrix inertia{};
            for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j) {
                const double identity=i==j?1:0;
                inertia[i][j]=tidal*(q[i]*q[j]-identity*dot(q,q)/3)
                    +spin*(omega[i]*omega[j]-identity*dot(omega,omega)/3
                        -n2*((i==2 && j==2?1:0)-identity/3));
            }
            // Park2021 eqs 30..34 and 54: inertia is normalized by m_M R_M^2.
            c[2][0]-=0.5*(2*inertia[2][2]-inertia[0][0]-inertia[1][1]);
            c[2][1]-=inertia[0][2];s[2][1]-=inertia[1][2];
            c[2][2]+=0.25*(inertia[1][1]-inertia[0][0]);s[2][2]-=0.5*inertia[0][1];
        }
        // Earth lies at -r in the Moon frame. Evaluate there explicitly:
        // odd degrees have different parity from the old quadrupole path.
        const Vector g=harmonic_acceleration<6>((-1.0)*body,options_.moon_radius_km,
            -(gm[0]+gm[1])*day2,options_.moon_degree,c,s);
        for(std::size_t j=0;j<3;++j)for(std::size_t i=0;i<3;++i)result[j]+=axes[i][j]*g[i];
        if(options_.external_figure_interactions)for(std::size_t source=2;source<=7;++source) {
            const auto q0=relative_positions[source]-r;
            const Vector point{dot(axes[0],q0),dot(axes[1],q0),dot(axes[2],q0)};
            const auto field=quadrupole_acceleration(point,options_.moon_radius_km,-gm[source]*day2,
                c[2][0],c[2][1],s[2][1],c[2][2],s[2][2]);
            for(std::size_t j=0;j<3;++j)for(std::size_t i=0;i<3;++i)result[j]+=axes[i][j]*field[i];
        }
    }
    return result;
}
Ephemeris::Ephemeris(spk::Kernel kernel,Instant<TDB> seed,Instant<TDB> first,Instant<TDB> last,
    IntegrationOptions integration,ForceOptions forces)
    :orbit_(Forces(kernel,forces).prepared(first,last),seed,seed_state(kernel,seed),first,last,integration){}
Ephemeris::Ephemeris(Forces forces,Instant<TDB> seed,State state,Instant<TDB> first,Instant<TDB> last,IntegrationOptions integration)
    :orbit_(forces.prepared(first,last),seed,state,first,last,integration){}
GeocentricEclipticJ2000State Ephemeris::operator()(Instant<TDB> t)const {
    return frames::mean_ecliptic_j2000_to_icrf().inverse()(orbit_(t));
}
GeocentricEclipticJ2000Position Ephemeris::position(Instant<TDB> t)const {
    return frames::mean_ecliptic_j2000_to_icrf().inverse()(orbit_.position(t));
}
} // namespace spp::model::lunar_numerical
