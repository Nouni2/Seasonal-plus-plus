// Numerical error/completeness gates coupled to Moon event-search timings.
#include "spp/spp.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "core/root_search.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <format>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
using spp::Time;
using spp::moon::HorizonKind;
using spp::model::elp_mpp02::CachedEphemeris;
using spp::model::elp_mpp02::Ephemeris;
using spp::model::elp_mpp02::Fit;
void require(bool condition,const char* why){if(!condition)throw std::runtime_error(why);}
double seconds(const Time& a,const Time& b) {
    const auto x=a.utc(),y=b.utc();return (static_cast<double>(x.day-y.day)+(x.frac-y.frac))*86400;
}
template<class F> double milliseconds(F f) {
    const auto a=Clock::now();f();return std::chrono::duration<double,std::milli>(Clock::now()-a).count();
}
template<class Exception,class F>void throws(F f){try{f();}catch(const Exception&){return;}throw std::runtime_error("expected search exception not raised");}

void scalar_cases() {
    auto run=[](const char* name,auto f,double curvature,int count,bool contact=false) {
        int calls=0;std::vector<spp::detail::ScalarRoot> roots;
        const auto cost=milliseconds([&]{roots=spp::detail::bounded_roots([&](double x){++calls;return f(x);},1,1e-8,1e-10,curvature);});
        if(static_cast<int>(roots.size())!=count) {
            std::fprintf(stderr,"%s: expected %d roots, found %zu (%d evaluations)\n",name,count,roots.size(),calls);
            for(const auto& r:roots)std::fprintf(stderr,"  x=%.17g f=%.17g rising=%d contact=%d\n",r.x,f(r.x),r.rising,r.contact);
            throw std::runtime_error("scalar root count wrong");
        }
        for(const auto& r:roots){require(r.contact==contact,"scalar contact classification wrong");require(std::abs(f(r.x))<2e-8,"scalar root residual too large");}
        std::printf("  %-25s %3zu roots  %5d evaluations  %.6f ms\n",name,roots.size(),calls,cost);
    };
    std::puts("Adaptive scalar search: analytic count/residual gates and measured work");
    run("two same-sign endpoints",[](double x){return (x-0.49)*(x-0.51);},2,2);
    run("subsecond window",[](double x){return (x-0.5)*(x-0.5)-1e-11;},2,2);
    run("tangent minimum",[](double x){return (x-0.5)*(x-0.5);},2,1,true);
    run("tangent maximum",[](double x){return -(x-0.5)*(x-0.5);},2,1,true);
    run("tangent near start",[](double x){return (x-0.0001)*(x-0.0001);},2,1,true);
    run("tangent near end",[](double x){return -(x-0.9999)*(x-0.9999);},2,1,true);
    run("contact at start",[](double x){return x*x;},2,1,true);
    run("contact at end",[](double x){return (x-1)*(x-1);},2,1,true);
    run("monotonic near start",[](double x){return 1e-5*(x+5e-6);},2,0);
    run("monotonic near end",[](double x){return 1e-5*(x-1-5e-6);},2,0);
    run("shallow noncontact",[](double x){return (x-0.5)*(x-0.5)+2e-10;},2,0);
    run("no contact",[](double x){return (x-0.5)*(x-0.5)+1e-5;},2,0);
    run("endpoint roots",[](double x){return x*(x-1);},2,2);
    throws<std::runtime_error>([]{(void)spp::detail::bounded_roots([](double){return 0.0;},1,1e-8,1e-10,2);});
    throws<std::invalid_argument>([]{(void)spp::detail::bounded_roots([](double x){return x;},1,0,1e-10,2);});
}

struct ReferenceEvent{HorizonKind kind;Time time;};
// Independent deliberately dense reference: a one-minute scan followed by
// ordinary bracket refinement. It is used for ordinary crossings, not as a
// completeness oracle for grazing or very short windows.
std::vector<ReferenceEvent> dense(const Time& start,const spp::Observer& obs,CachedEphemeris& model) {
    const spp::Atmosphere airless{.refraction=false};
    const auto atmosphere=spp::sky::at_elevation({},obs.elev_m);
    auto place=[&](double x){return spp::moon::coordinates(start.shifted(x),obs,model,airless);};
    auto alt=[&](double x){const auto c=place(x);return spp::sky::refraction(c.geometric_alt_deg+c.semidiameter_deg,atmosphere);};
    auto meridian=[&](double x){return std::sin(place(x).hour_angle_deg*spp::deg2rad);};
    std::vector<ReferenceEvent> result;
    double previous=0;
    auto c=place(0);
    double fa=spp::sky::refraction(c.geometric_alt_deg+c.semidiameter_deg,atmosphere),ha=std::sin(c.hour_angle_deg*spp::deg2rad);
    for(int i=1;i<=1440;++i) {
        const double x=static_cast<double>(i)/1440;c=place(x);
        const double fb=spp::sky::refraction(c.geometric_alt_deg+c.semidiameter_deg,atmosphere),hb=std::sin(c.hour_angle_deg*spp::deg2rad);
        auto bracket=[&](auto function,double a,double b) {
            if(a==0)return previous;
            if(b==0)return x;
            return a>0?spp::detail::refine_root(function,previous,a,x,b,0.002/86400,100)
                :spp::detail::refine_root(function,x,b,previous,a,0.002/86400,100);
        };
        if(fa==0 || (fa<0 && fb>0)||(fa>0 && fb<0)) {
            const double root=bracket(alt,fa,fb);
            if(root<1)result.push_back({fb>fa?HorizonKind::Rise:HorizonKind::Set,start.shifted(root)});
        }
        if(ha==0 || (ha<0 && hb>0)||(ha>0 && hb<0)) {
            const double root=bracket(meridian,ha,hb);
            if(root<1)result.push_back({std::cos(place(root).hour_angle_deg*spp::deg2rad)>0?
                HorizonKind::UpperTransit:HorizonKind::LowerTransit,start.shifted(root)});
        }
        previous=x;fa=fb;ha=hb;
    }
    std::sort(result.begin(),result.end(),[](const auto& a,const auto& b){return seconds(a.time,b.time)<0;});
    return result;
}

void comparison() {
    struct Case {int month,day; spp::Observer obs;const char* name;};
    const std::vector<Case> cases={{11,3,{48.8566,2.3522,35},"Paris"},{11,3,{-33.8688,151.2093,58},"Sydney"},
        {3,20,{0,0,0},"equator"},{6,21,{78,15,0},"Arctic"},{6,21,{-78,15,0},"Antarctic"},
        {11,3,{90,0,0},"north pole"},{11,3,{-90,0,0},"south pole"}};
    double max_time=0,max_altitude=0,max_meridian=0;
    for(Fit fit:{Fit::DE405,Fit::LLR}) {
        std::printf("Moon event search, %s (cold models, one UTC day; then warm repeat):\n",fit==Fit::DE405?"DE405":"LLR");
        for(const auto& q:cases) {
            const auto start=Time::from_utc(2026,q.month,q.day,0,0,0,0);
            CachedEphemeris adaptive(fit),reference(fit);
            spp::moon::HorizonResult found;
            const double cold=milliseconds([&]{found=spp::moon::horizon(start,start.shifted(1),q.obs,adaptive);});
            const auto calls=adaptive.source_evaluations();
            std::vector<ReferenceEvent> expected;
            const double scan=milliseconds([&]{expected=dense(start,q.obs,reference);});
            require(found.events.size()==expected.size(),"adaptive/dense event count differs");
            require(found.duration_hours>=0 && found.duration_hours<=24.000001,"above-horizon duration out of range");
            for(std::size_t i=0;i<found.events.size();++i) {
                const auto& e=found.events[i];require(e.kind==expected[i].kind,"adaptive/dense event ordering differs");
                require(seconds(e.time,start)>=0 && seconds(e.time,start.shifted(1))<0,"event outside half-open interval");
                max_time=std::max(max_time,std::abs(seconds(e.time,expected[i].time)));
                const auto c=spp::moon::coordinates(e.time,q.obs,adaptive,{.refraction=false});
                if(e.kind==HorizonKind::Rise || e.kind==HorizonKind::Set) {
                    const auto atm=spp::sky::at_elevation({},q.obs.elev_m);
                    max_altitude=std::max(max_altitude,std::abs(spp::sky::refraction(c.geometric_alt_deg+c.semidiameter_deg,atm))*3600);
                } else max_meridian=std::max(max_meridian,std::abs(std::sin(c.hour_angle_deg*spp::deg2rad)));
            }
            std::vector<double> repeat;
            for(int i=0;i<3;++i)repeat.push_back(milliseconds([&]{auto result=spp::moon::horizon(start,start.shifted(1),q.obs,adaptive);require(result.events.size()==found.events.size(),"repeat changed event count");}));
            std::sort(repeat.begin(),repeat.end());
            std::printf("  %-12s %zu events  above %7.4f h  adaptive %8.3f ms  dense %8.3f ms  warm %7.3f ms  %llu ELP states\n",
                q.name,found.events.size(),found.duration_hours,cold,scan,repeat[1],static_cast<unsigned long long>(calls));
            if(fit==Fit::DE405 && q.name==std::string_view("Paris")) {
                const Ephemeris direct(fit);spp::moon::HorizonResult exact;
                const double direct_cost=milliseconds([&]{exact=spp::moon::horizon(start,start.shifted(1),q.obs,direct);});
                require(exact.events.size()==found.events.size(),"direct/cached event count differs");
                for(std::size_t i=0;i<exact.events.size();++i)
                    require(exact.events[i].kind==found.events[i].kind && std::abs(seconds(exact.events[i].time,found.events[i].time))<0.05,
                            "interpolation changed horizon events");
                std::printf("    direct full-series search %.3f ms; cached/direct event times agree within 0.05 s\n",direct_cost);
            }
        }
    }
    std::printf("  max adaptive/dense difference %.9f s; limb altitude residual %.9f arcsec; sin(hour angle) residual %.9g\n",
        max_time,max_altitude,max_meridian);
    require(max_time<0.05 && max_altitude<0.5 && max_meridian<1e-6,"event numerical error gate failed");
}

void grazing_and_boundaries() {
    const auto start=Time::from_utc(2026,11,3,0,0,0,0);
    const spp::Observer obs{48.8566,2.3522,35};CachedEphemeris model;
    const spp::Atmosphere airless{.refraction=false};
    auto alt=[&](double d){return spp::moon::horizontal(start.shifted(d),obs,model,airless).alt_deg;};
    int peak=1;double highest=alt(1.0/24);
    for(int i=2;i<24;++i)if(const double value=alt(i/24.0);value>highest){peak=i;highest=value;}
    require(peak>0 && peak<23,"grazing workload has no interior maximum");
    double a=(peak-1)/24.0,b=(peak+1)/24.0;
    constexpr double ratio=0.6180339887498948482;
    double x=b-ratio*(b-a),y=a+ratio*(b-a),fx=alt(x),fy=alt(y);
    while(b-a>1e-10) {
        if(fx>fy){b=y;y=x;fy=fx;x=b-ratio*(b-a);fx=alt(x);}
        else{a=x;x=y;fx=fy;y=a+ratio*(b-a);fy=alt(y);}
    }
    const double maximum=std::max(fx,fy);
    const spp::moon::HorizonOptions touch{.altitude_deg=maximum,.upper_limb=false};
    auto contacts=spp::moon::horizon(start,start.shifted(1),obs,model,airless,touch);
    int grazes=0,crossings=0;
    for(const auto& e:contacts.events) {
        grazes+=e.kind==HorizonKind::Graze;
        crossings+=e.kind==HorizonKind::Rise || e.kind==HorizonKind::Set;
    }
    require(grazes==1 && crossings==0,"stationary horizon contact misclassified");
    auto short_window=touch;short_window.altitude_deg-=1e-5;
    const auto pair=spp::moon::horizon(start,start.shifted(1),obs,model,airless,short_window);
    crossings=0;for(const auto& e:pair.events)crossings+=e.kind==HorizonKind::Rise || e.kind==HorizonKind::Set;
    require(crossings==2 && pair.duration_hours>0 && pair.duration_hours<1.0/60,"short same-sign window was missed");
    auto missed=touch;missed.altitude_deg+=1e-5;
    const auto below=spp::moon::horizon(start,start.shifted(1),obs,model,airless,missed);
    for(const auto& e:below.events)require(e.kind!=HorizonKind::Graze && e.kind!=HorizonKind::Rise && e.kind!=HorizonKind::Set,"spurious horizon contact");

    const auto top_time=start.shifted(fx>fy?x:y);
    std::puts("Boundary-contact visibility/duration gates:");
    const auto after_max=spp::moon::horizon(top_time,top_time.shifted(1.0/24),obs,model,airless,touch);
    require(!after_max.above_at_start && after_max.duration_hours<1e-6,"boundary maximum changed visibility");
    int trough=1;double lowest=alt(1.0/24);
    for(int i=2;i<24;++i)if(const double value=alt(i/24.0);value<lowest){trough=i;lowest=value;}
    require(trough>0 && trough<23,"boundary workload has no interior minimum");
    a=(trough-1)/24.0;b=(trough+1)/24.0;
    x=b-ratio*(b-a);y=a+ratio*(b-a);fx=alt(x);fy=alt(y);
    while(b-a>1e-10) {
        if(fx<fy){b=y;y=x;fy=fx;x=b-ratio*(b-a);fx=alt(x);}
        else{a=x;x=y;fx=fy;y=a+ratio*(b-a);fy=alt(y);}
    }
    const auto bottom_time=start.shifted(fx<fy?x:y);
    const auto after_min=spp::moon::horizon(bottom_time,bottom_time.shifted(1.0/24),obs,model,airless,
        {.altitude_deg=std::min(fx,fy),.upper_limb=false});
    require(after_min.above_at_start && std::abs(after_min.duration_hours-1)<1e-6,"boundary minimum lost its above-horizon duration");

    const auto moment=Time::from_utc(2026,11,3,12,0,0,0);
    const double target=spp::moon::horizontal(moment,obs,model,airless).alt_deg;
    const spp::moon::HorizonOptions edge{.altitude_deg=target,.upper_limb=false};
    const auto right=spp::moon::horizon(moment,moment.shifted(0.125),obs,model,airless,edge);
    bool at_start=false;
    for(const auto& e:right.events)if((e.kind==HorizonKind::Rise || e.kind==HorizonKind::Set) && std::abs(seconds(e.time,moment))<0.01)at_start=true;
    require(at_start,"inclusive start missed its crossing");
    const auto left=spp::moon::horizon(moment.shifted(-0.125),moment,obs,model,airless,edge);
    for(const auto& e:left.events)require(seconds(e.time,moment)<0,"exclusive end included its crossing");
    require(spp::moon::horizon(moment,moment,obs,model).events.empty(),"empty interval not empty");
    throws<std::invalid_argument>([&]{spp::moon::horizon(moment,moment.shifted(-1),obs,model);});
    throws<std::invalid_argument>([&]{spp::moon::horizon(moment,moment.shifted(1),obs,model,{}, {.time_tolerance_s=0});});
    std::printf("Grazing/short-window/boundary gates: PASS; narrow visibility %.6f seconds\n",pair.duration_hours*3600);
}
}  // namespace

int main() {
    try {scalar_cases();comparison();grazing_and_boundaries();std::puts("PASS: Moon event performance and numerical gates");return 0;}
    catch(const std::exception& e){std::fprintf(stderr,"FAIL Moon event benchmark: %s\n",e.what());return 1;}
}
