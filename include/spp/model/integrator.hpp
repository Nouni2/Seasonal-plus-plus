#pragma once

#include "spp/core/state.hpp"
#include "spp/core/time.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace spp::model {

struct IntegrationOptions {
    double relative_tolerance=1e-12;
    double position_tolerance_km=1e-7;
    double velocity_tolerance_km_day=1e-7;
    double initial_step_days=1.0/144;
    double maximum_step_days=1.0/48;
    double minimum_step_days=1e-8;
    std::size_t maximum_steps=250000;
};
struct IntegrationStats {std::size_t accepted_steps=0,rejected_steps=0,force_evaluations=0;};

// Explicit finite-window Dormand-Prince 5(4) propagation. Force(t,state)
// returns acceleration in km/TDB day^2. Origin/frame are retained in the type.
// Integrates from one seed in both directions; never reseeds from an ephemeris.
// Const queries use binary lookup and quintic Hermite dense position output;
// returned velocity is its analytic derivative. No force calls on queries.
template<class Force,class Origin=origin::Geocentric,class Frame=frame::Icrf>
class IntegratedOrbit {
public:
    using State=CartesianState<Origin,Frame>;
    using Position=CartesianPosition<Origin,Frame>;
    using Epoch=Instant<TDB>;
    IntegratedOrbit(Force force,Epoch seed,State state,Epoch first,Epoch last,IntegrationOptions options={})
        :force_(std::move(force)),seed_(normalize(seed)),initial_(state),options_(options) {
        first_=difference(normalize(first),seed_);last_=difference(normalize(last),seed_);
        validate();
        if(first_<0)propagate(first_);
        if(last_>0)propagate(last_);
        std::sort(segments_.begin(),segments_.end(),[](const auto& a,const auto& b){return a.low<b.low;});
    }
    State operator()(Epoch t)const {
        const double x=difference(normalize(t),seed_);
        if(x<first_ || x>last_)throw std::out_of_range("spp integrator: query outside propagated interval");
        if(x==0)return initial_;
        auto i=std::upper_bound(segments_.begin(),segments_.end(),x,[](double value,const Segment& s){return value<s.low;});
        if(i==segments_.begin())throw std::out_of_range("spp integrator: no segment");
        --i;const double u=(x-i->start)/i->step;
        State result{};
        for(std::size_t axis=0;axis<3;++axis) {
            const auto& c=i->coefficients[axis];
            result.position_km[axis]=c[0]+u*(c[1]+u*(c[2]+u*(c[3]+u*(c[4]+u*c[5]))));
            result.velocity_km_per_day[axis]=(c[1]+u*(2*c[2]+u*(3*c[3]+u*(4*c[4]+u*5*c[5]))))/i->step;
        }
        return result;
    }
    Position position(Epoch t)const{return {(*this)(t).position_km};}
    const IntegrationStats& stats()const{return stats_;}
    const Force& force()const{return force_;}
    std::size_t storage_bytes()const{return segments_.capacity()*sizeof(Segment);}
private:
    using Y=std::array<double,6>;
    struct Segment {double low,start,step;std::array<std::array<double,6>,3> coefficients;};
    Force force_;Epoch seed_;State initial_;IntegrationOptions options_;IntegrationStats stats_;
    double first_=0,last_=0;std::vector<Segment> segments_;
    static Epoch normalize(Epoch t) {
        if(!std::isfinite(t.frac))throw std::invalid_argument("spp integrator: non-finite TDB date");
        if(t.day<-1000000000 || t.day>1000000000 || std::abs(t.frac)>1000000000)
            throw std::out_of_range("spp integrator: unusable TDB date");
        const double whole=std::floor(t.frac);return {t.day+static_cast<std::int64_t>(whole),t.frac-whole};
    }
    static double difference(Epoch a,Epoch b){return static_cast<double>(a.day-b.day)+(a.frac-b.frac);}
    void validate()const {
        const auto& o=options_;
        for(double x:{o.relative_tolerance,o.position_tolerance_km,o.velocity_tolerance_km_day,
            o.initial_step_days,o.maximum_step_days,o.minimum_step_days})
            if(!std::isfinite(x) || x<=0)throw std::invalid_argument("spp integrator: positive finite tolerances/steps required");
        if(o.relative_tolerance<4*std::numeric_limits<double>::epsilon() || o.relative_tolerance>0.01
            || o.minimum_step_days>o.maximum_step_days || o.initial_step_days<o.minimum_step_days
            || o.initial_step_days>o.maximum_step_days || !o.maximum_steps || o.maximum_steps>10000000
            || first_>0 || last_<0 || first_>last_ || last_-first_>1000000)
            throw std::invalid_argument("spp integrator: invalid limits or seed outside interval");
        for(std::size_t i=0;i<3;++i)if(!std::isfinite(initial_.position_km[i]) || !std::isfinite(initial_.velocity_km_per_day[i]))
            throw std::invalid_argument("spp integrator: non-finite seed state");
    }
    Y derivative(double x,const Y& y) {
        State state{{y[0],y[1],y[2]},{y[3],y[4],y[5]}};
        for(double v:y)if(!std::isfinite(v))throw std::runtime_error("spp integrator: non-finite trial state");
        const auto a=force_(normalize({seed_.day,seed_.frac+x}),state);++stats_.force_evaluations;
        Y result{y[3],y[4],y[5],a[0],a[1],a[2]};
        for(double v:result)if(!std::isfinite(v))throw std::runtime_error("spp integrator: non-finite force");
        return result;
    }
    void store(double x,double h,const Y& y,const Y& end,const Y& f,const Y& fend) {
        Segment s{std::min(x,x+h),x,h,{}};
        for(std::size_t j=0;j<3;++j) {
            auto& c=s.coefficients[j];c[0]=y[j];c[1]=h*y[j+3];c[2]=0.5*h*h*f[j+3];
            const double p=end[j]-y[j]-c[1]-c[2],v=h*end[j+3]-c[1]-2*c[2],a=h*h*fend[j+3]-2*c[2];
            c[3]=10*p-4*v+0.5*a;c[4]=-15*p+7*v-a;c[5]=6*p-3*v+0.5*a;
        }
        segments_.push_back(s);
    }
    void propagate(double limit) {
        // Dormand & Prince (1980), embedded fifth/fourth-order tableau.
        constexpr double a[7][7]={{0},{1.0/5},{3.0/40,9.0/40},{44.0/45,-56.0/15,32.0/9},
            {19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729},
            {9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656},
            {35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84}};
        constexpr double c[]={0,1.0/5,3.0/10,4.0/5,8.0/9,1,1};
        constexpr double e[]={35.0/384-5179.0/57600,0,500.0/1113-7571.0/16695,
            125.0/192-393.0/640,-2187.0/6784+92097.0/339200,11.0/84-187.0/2100,-1.0/40};
        Y y{initial_.position_km[0],initial_.position_km[1],initial_.position_km[2],
            initial_.velocity_km_per_day[0],initial_.velocity_km_per_day[1],initial_.velocity_km_per_day[2]};
        double x=0,h=std::copysign(options_.initial_step_days,limit);std::array<Y,7> k{};
        k[0]=derivative(x,y);bool rejected=false;
        while(x!=limit) {
            if(stats_.accepted_steps+stats_.rejected_steps>=options_.maximum_steps)
                throw std::runtime_error("spp integrator: step budget exhausted");
            const double remaining=limit-x;
            if(std::abs(h)>std::abs(remaining))h=remaining;
            if(x+h==x || (std::abs(h)<options_.minimum_step_days && std::abs(remaining)>options_.minimum_step_days))
                throw std::runtime_error("spp integrator: step underflow before tolerance reached");
            Y end{};
            for(std::size_t s=1;s<7;++s) {
                Y stage=y;
                for(std::size_t j=0;j<6;++j) {
                    double sum=0;for(std::size_t q=0;q<s;++q)sum+=a[s][q]*k[q][j];
                    stage[j]+=h*sum;
                }
                k[s]=derivative(s==6?x+h:x+c[s]*h,stage);
                if(s==6)end=stage;
            }
            double error=0;
            for(std::size_t j=0;j<6;++j) {
                double residual=0;for(std::size_t s=0;s<7;++s)residual+=e[s]*k[s][j];
                const double absolute=j<3?options_.position_tolerance_km:options_.velocity_tolerance_km_day;
                const double scale=absolute+options_.relative_tolerance*std::max(std::abs(y[j]),std::abs(end[j]));
                error=std::max(error,std::abs(h*residual)/scale);
            }
            if(!std::isfinite(error))throw std::runtime_error("spp integrator: non-finite error estimate");
            double factor=error==0?5:std::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
            if(error<=1) {
                store(x,h,y,end,k[0],k[6]);++stats_.accepted_steps;
                x=(h==remaining)?limit:x+h;y=end;k[0]=k[6];
                if(rejected)factor=std::min(1.0,factor);
                rejected=false;
            } else {++stats_.rejected_steps;rejected=true;}
            h=std::copysign(std::min(options_.maximum_step_days,std::abs(h)*factor),limit);
        }
    }
};
} // namespace spp::model
