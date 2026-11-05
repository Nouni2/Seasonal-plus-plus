// Lunar phase and observer application; OpenGL 3.3, no network/runtime assets.
#include "gl_view.hpp"
#include <GLFW/glfw3.h>
#include "offscreen.hpp"
#include "task_worker.hpp"
#include "spp/spp.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "spp/model/lunar_numerical.hpp"
#include "spp/model/meeus.hpp"

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <format>
#include <memory>
#include <mutex>
#include <optional>
#include <print>
#include <string>
#include <thread>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
using spp::Time;
using spp::model::elp_mpp02::Fit;
using spp::app::Rect;
using spp::app::ink;
using spp::app::muted;
using spp::app::accent;

Time now() {
    const auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    return Time::from_jd(2440587,0.5+static_cast<double>(ms)/86400000.0);
}
double days(const Time& a,const Time& b) {
    const auto x=a.utc(),y=b.utc();return static_cast<double>(x.day-y.day)+(x.frac-y.frac);
}
std::string stamp(const Time& t) {
    auto s=t.iso_utc();std::replace(s.begin(),s.end(),'T',' ');if(s.back()=='Z')s.pop_back();return s;
}
spp::CivilDateTime civil_timestamp(const std::string& s) {
    int y=0,m=0,d=0,h=0,mi=0,n=0;double sec=0;
    if(std::sscanf(s.c_str(),"%d-%d-%d %d:%d:%lf%n",&y,&m,&d,&h,&mi,&sec,&n)!=6
        || n!=static_cast<int>(s.size()) || y< -2000 || y>4000 || !spp::is_valid({y,m,d})
        || h<0 || h>23 || mi<0 || mi>59 || !std::isfinite(sec) || sec<0 || sec>=60)
        throw std::invalid_argument("Use a valid YYYY-MM-DD HH:MM:SS timestamp.");
    return {y,m,d,h,mi,sec};
}
Time timestamp(const std::string& s){return Time::from_civil(civil_timestamp(s));}
using Numerical=spp::model::lunar_numerical::Ephemeris;
double number(std::string_view s) {
    double x=0;const auto [end,error]=std::from_chars(s.data(),s.data()+s.size(),x);
    if(error!=std::errc{} || end!=s.data()+s.size() || !std::isfinite(x))
        throw std::invalid_argument("Enter a finite number.");
    return x;
}
struct Options {
    Time time=now();spp::Observer observer{48.8566,2.3522,35};Fit fit=Fit::DE405;
    int frames=0;bool offscreen=false,direct=false,refraction=true,events=false,grid=false;
    double step_seconds=60;std::string capture;
    std::shared_ptr<const Numerical> numerical;
};
Options options(int argc,char** argv) {
    Options o;std::vector<std::string> positional;
    bool integrated=false,fit_given=false;std::string kernel;std::optional<spp::Instant<spp::TDB>> seed;double span=90;
    for(int i=1;i<argc;++i) {
        const std::string s=argv[i];
        if(s=="--help") {
            std::println("usage: spp_moon [YYYY-MM-DD HH:MM:SS [LAT LON [ELEV]]] [--fit=de405|llr]\n"
                "  --direct  --no-refraction  --events  --grid\n"
                "  --model=integrated --kernel=PATH --seed=YYYY-MM-DDTHH:MM[:SS] (TDB), --span-days=90\n"
                "  --benchmark=N [--offscreen] [--step-seconds=60] [--capture=FILE.ppm]\n"
                "Keys: Space play/pause, arrows hour (+Shift day), +/- speed, N now, F fit,\n"
                "      O north/zenith up, D direct/cache, R refraction, E earthshine, G lit area, H daily sky, L libration grid.\n"
                "Click UTC/location fields to edit; Enter accepts, Escape cancels. Drag the timeline.");
            std::exit(0);
        } else if(s=="--offscreen")o.offscreen=true;
        else if(s=="--events")o.events=true;
        else if(s=="--grid")o.grid=true;
        else if(s=="--direct")o.direct=true;
        else if(s=="--no-refraction")o.refraction=false;
        else if(s=="--fit=llr"){o.fit=Fit::LLR;fit_given=true;}
        else if(s=="--fit=de405"){o.fit=Fit::DE405;fit_given=true;}
        else if(s=="--model=elp")integrated=false;
        else if(s=="--model=integrated")integrated=true;
        else if(s.starts_with("--kernel="))kernel=s.substr(9);
        else if(s.starts_with("--span-days="))span=number(std::string_view(s).substr(12));
        else if(s.starts_with("--seed=")) {
            std::string clock=s.substr(7);std::replace(clock.begin(),clock.end(),'T',' ');
            if(std::count(clock.begin(),clock.end(),':')==1)clock+=":00";
            const auto jd=spp::model::meeus::jd_from_civil(civil_timestamp(clock));seed={{jd.day,jd.frac}};
        }
        else if(s.starts_with("--benchmark=")) {
            const auto v=number(std::string_view(s).substr(12));
            if(v<1 || v>100000 || v!=std::floor(v))throw std::invalid_argument("benchmark frames must be 1..100000");
            o.frames=static_cast<int>(v);
        } else if(s.starts_with("--step-seconds="))o.step_seconds=number(std::string_view(s).substr(15));
        else if(s.starts_with("--capture="))o.capture=s.substr(10);
        else if(s.starts_with("--"))throw std::invalid_argument("unknown option: "+s);
        else positional.push_back(s);
    }
    if(!positional.empty()) {
        if(positional.size()!=2 && positional.size()!=4 && positional.size()!=5)
            throw std::invalid_argument("provide date/time and, optionally, latitude longitude elevation");
        std::string clock=positional[1];if(std::count(clock.begin(),clock.end(),':')==1)clock+=":00";
        o.time=timestamp(positional[0]+" "+clock);
        if(positional.size()>=4)o.observer={number(positional[2]),number(positional[3]),positional.size()==5?number(positional[4]):0};
    }
    if(o.offscreen && !o.frames)throw std::invalid_argument("--offscreen requires an explicit --benchmark=N");
    if(std::abs(o.step_seconds)>86400)throw std::invalid_argument("benchmark step must be within +/-86400 seconds");
    if(span<=0 || span>365)throw std::invalid_argument("--span-days must be positive and at most 365");
    if(integrated) {
        if(kernel.empty() || !seed)throw std::invalid_argument("integrated Moon requires --kernel and explicit TDB --seed");
        if(fit_given || o.direct)throw std::invalid_argument("--fit and --direct select ELP and cannot be combined with the integrated model");
        const auto start=Clock::now();
        o.numerical=std::make_shared<Numerical>(spp::model::spk::Kernel(kernel),*seed,
            spp::Instant<spp::TDB>{seed->day,seed->frac-span},spp::Instant<spp::TDB>{seed->day,seed->frac+span});
        const auto& stats=o.numerical->stats();
        std::println("Numerical Moon: {:.3f} ms propagation, {} accepted / {} rejected steps, {} force evaluations; +/-{:g} days from TDB seed.",
            std::chrono::duration<double,std::milli>(Clock::now()-start).count(),stats.accepted_steps,stats.rejected_steps,stats.force_evaluations,span);
        std::println("  Seed JD {} + {:.12f} TDB; prescribed EMB/planets, 1PN, Earth degree 5, lunar degree 6, tidal deformation and solar effects.",seed->day,seed->frac);
    } else if(!kernel.empty() || seed || span!=90)throw std::invalid_argument("kernel/seed/span options require --model=integrated");
    return o;
}

void draw_grid(spp::app::View& view,float cx,float cy,float radius,const spp::moon::Surface& surface,double angle) {
    using V=std::array<double,3>;
    struct Edge {V a,b;bool major;};
    static const auto mesh=[] {
        std::vector<Edge> edges;edges.reserve(1980);
        auto point=[](int longitude,int latitude) {
            const double l=longitude*spp::deg2rad,b=latitude*spp::deg2rad;
            return V{std::cos(b)*std::cos(l),std::cos(b)*std::sin(l),std::sin(b)};
        };
        for(int lat:{-60,-30,0,30,60})for(int lon=0;lon<360;lon+=2)
            edges.push_back({point(lon,lat),point(lon+2,lat),lat==0});
        for(int lon=0;lon<360;lon+=30)for(int lat=-90;lat<90;lat+=2)
            edges.push_back({point(lon,lat),point(lon,lat+2),lon==0});
        return edges;
    }();
    const double c=std::cos(angle),s=std::sin(angle);
    auto project=[&](V p) {
        V q{};
        for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)q[i]+=surface.body_to_view[i][j]*p[j];
        return V{c*q[0]+s*q[1],-s*q[0]+c*q[1],q[2]};
    };
    for(const auto& edge:mesh) {
        V a=project(edge.a),b=project(edge.b);
        if(a[2]<=0 && b[2]<=0)continue;
        if(a[2]<0 || b[2]<0) {
            const double u=a[2]/(a[2]-b[2]);
            V limb{a[0]+u*(b[0]-a[0]),a[1]+u*(b[1]-a[1]),0};
            const double norm=std::hypot(limb[0],limb[1]);
            if(norm>0){limb[0]/=norm;limb[1]/=norm;}
            if(a[2]<0)a=limb;else b=limb;
        }
        const auto color=edge.major?spp::app::Color{0.82f,0.69f,0.43f,0.65f}:spp::app::Color{0.48f,0.55f,0.64f,0.48f};
        view.line(cx+radius*static_cast<float>(a[0]),cy-radius*static_cast<float>(a[1]),
            cx+radius*static_cast<float>(b[0]),cy-radius*static_cast<float>(b[1]),color,edge.major?1.2f:0.8f);
    }
    const auto origin=project({1,0,0});
    if(origin[2]>0) {
        const float x=cx+radius*static_cast<float>(origin[0]),y=cy-radius*static_cast<float>(origin[1]);
        view.line(x-4,y,x+4,y,accent,1.2f);view.line(x,y-4,x,y+4,accent,1.2f);
        view.text(x+6,y-13,"ME 0,0",0.58f,accent);
    }
}

struct Calendar {
    Time anchor;Fit fit;std::uint64_t generation;
    std::vector<double> fraction;
    std::vector<spp::moon::PhaseEvent> events;
    std::string error;
    std::shared_ptr<const Numerical> numerical;
};
using Cancelled=spp::app::WorkerCancelled;
template<class Result>using Worker=spp::app::TaskWorker<Result>;
template<class Cancel>
struct CancellableModel {
    spp::model::elp_mpp02::CachedEphemeris source;
    Cancel cancelled;
    std::shared_ptr<const Numerical> numerical;
    spp::moon::Position position(spp::Instant<spp::TDB> t) {
        if(cancelled())throw Cancelled{};
        return numerical?numerical->position(t):source.position(t);
    }
};
template<class Cancel>
void fill_calendar(Calendar& out,Cancel cancelled) {
    CancellableModel<Cancel> model{spp::model::elp_mpp02::CachedEphemeris(out.fit),cancelled,out.numerical};
    out.events=spp::moon::phases(out.anchor.shifted(-32),out.anchor.shifted(74),model);
    for(int i=0;i<=112;++i) {
        if(cancelled())throw Cancelled{};
        out.fraction.push_back(spp::moon::phase(out.anchor.shifted(i*42.0/112),model).illuminated_fraction);
    }
}

struct DailySky {
    Time anchor;
    spp::Observer observer;
    Fit fit;
    bool refraction;
    std::uint64_t generation=0;
    spp::moon::HorizonResult horizon;
    std::vector<double> altitude;
    std::string error;
    std::shared_ptr<const Numerical> numerical;
};
template<class Cancel>
void fill_calendar(DailySky& out,Cancel cancelled) {
    CancellableModel<Cancel> model{spp::model::elp_mpp02::CachedEphemeris(out.fit),cancelled,out.numerical};
    const spp::Atmosphere atm{.refraction=out.refraction};
    out.horizon=spp::moon::horizon(out.anchor,out.anchor.shifted(1),out.observer,model,atm);
    for(int i=0;i<=48;++i)
        out.altitude.push_back(spp::moon::horizontal(out.anchor.shifted(i/48.0),out.observer,model,atm).alt_deg);
}

class Application {
    struct Direct {
        spp::model::elp_mpp02::Ephemeris model;std::uint64_t calls=0;
        explicit Direct(Fit fit):model(fit){}
        spp::moon::Position position(spp::Instant<spp::TDB> t){++calls;return model.position(t);}
    };
    spp::app::View view_;
    spp::model::elp_mpp02::CachedEphemeris cached_;
    Direct direct_;
    std::unique_ptr<Worker<Calendar>> worker_;
    std::unique_ptr<Worker<DailySky>> daily_worker_;
    std::shared_ptr<const Calendar> calendar_;
    std::optional<Time> requested_anchor_;std::uint64_t request_id_=0;
    std::shared_ptr<const DailySky> daily_;
    std::optional<DailySky> requested_day_;
    std::uint64_t daily_id_=0,calendar_updates_=0,daily_updates_=0;
    spp::moon::Observation observation_{};
    bool dirty_=true,valid_=false,pressed_=false,down_=false,select_all_=false;
    bool north_up_=true,earthshine_=false,area_=false;
    int edit_=-1;std::string buffer_,error_;
    double mx_=0,my_=0;
    float width_=1280,height_=900;
    Rect timeline_{},day_chart_{};
    double calculation_us_=0;

    bool button(Rect rect,std::string_view label,bool active=false) {
        const bool hovered=rect.contains(mx_,my_);
        view_.quad(rect,active?spp::app::Color{0.22f,0.19f,0.13f,1}:
            hovered?spp::app::Color{0.15f,0.18f,0.22f,1}:spp::app::Color{0.085f,0.10f,0.13f,1});
        view_.text(rect.x+8,rect.y+3,label,0.78f,active?accent:ink);
        if(pressed_ && hovered){pressed_=false;edit_=-1;if(worker_)glfwPostEmptyEvent();return true;}return false;
    }
    void field(int id,Rect rect,std::string text) {
        const bool focus=edit_==id;
        view_.quad(rect,focus?spp::app::Color{0.16f,0.18f,0.22f,1}:spp::app::Color{0.07f,0.08f,0.11f,1});
        view_.text(rect.x+5,rect.y+3,focus?buffer_:text,0.78f,focus?accent:ink);
        if(pressed_ && rect.contains(mx_,my_)) {
            edit_=id;buffer_=std::move(text);select_all_=true;pressed_=false;
            if(worker_)glfwPostEmptyEvent();
        }
    }
    void change_fit() {
        if(config.numerical)return;
        config.fit=config.fit==Fit::DE405?Fit::LLR:Fit::DE405;
        cached_=spp::model::elp_mpp02::CachedEphemeris(config.fit);
        direct_=Direct(config.fit);
        calendar_.reset();requested_anchor_.reset();dirty_=true;
    }
    void update_calendar() {
        if(worker_) {
            const auto result=worker_->ready();
            if(result && result->generation==request_id_)calendar_=result;
        }
        if(!requested_anchor_ || days(config.time,*requested_anchor_)<3 || days(config.time,*requested_anchor_)>39) {
            const auto u=config.time.utc();
            const auto midnight=Time::from_jd(u.day-(u.frac<0.5?1:0),0.5);
            requested_anchor_=midnight.shifted(-21);
            ++calendar_updates_;
            calendar_.reset();
            if(worker_)request_id_=worker_->request(Calendar{*requested_anchor_,config.fit,0,{},{},{},config.numerical});
            else {
                auto result=std::make_shared<Calendar>(Calendar{*requested_anchor_,config.fit,0,{},{},{},config.numerical});
                try {fill_calendar(*result,[]{return false;});}catch(const std::exception& e){result->error=e.what();}
                calendar_=std::move(result);
            }
        }
    }
    void update_daily() {
        if(!config.events)return;
        if(daily_worker_) {
            const auto result=daily_worker_->ready();
            if(result && result->generation==daily_id_)daily_=result;
        }
        const auto u=config.time.utc();
        const auto anchor=Time::from_jd(u.day-(u.frac<0.5?1:0),0.5);
        const auto& o=config.observer;
        if(!requested_day_ || days(anchor,requested_day_->anchor)!=0
            || o.lat_deg!=requested_day_->observer.lat_deg || o.lon_deg!=requested_day_->observer.lon_deg
            || o.elev_m!=requested_day_->observer.elev_m || config.fit!=requested_day_->fit
            || config.refraction!=requested_day_->refraction) {
            requested_day_=DailySky{anchor,o,config.fit,config.refraction,0,{},{},{},config.numerical};
            ++daily_updates_;
            daily_.reset();
            if(daily_worker_)daily_id_=daily_worker_->request(*requested_day_);
            else {
                auto result=std::make_shared<DailySky>(*requested_day_);
                try {fill_calendar(*result,[]{return false;});}catch(const std::exception& e){result->error=e.what();}
                daily_=std::move(result);
            }
        }
    }
    void daily_panel(float top) {
        if(!daily_ || !daily_->error.empty()) {
            view_.text(24,top+55,daily_?daily_->error:"Computing daily sky...",0.78f,muted);
            return;
        }
        const auto& h=daily_->horizon;
        view_.text(395,top+2,std::format("{} UTC | Above horizon {:.3f} h",
            spp::iso_date(daily_->anchor.to_civil().date()),h.duration_hours),0.73f,muted);
        view_.text(24,top+29,"Center altitude; rise/set use the upper limb. Click an event or drag the chart.",0.62f,muted);
        day_chart_={54,top+52,width_-100,60};
        const float bottom=day_chart_.y+day_chart_.h;
        auto level=[&](double altitude){return day_chart_.y+day_chart_.h*static_cast<float>((90-altitude)/180);};
        for(int altitude:{-60,0,60}) {
            const float y=level(altitude);
            view_.text(24,y-7,std::format("{:+d}",altitude),0.52f,muted);
            view_.line(day_chart_.x,y,day_chart_.x+day_chart_.w,y,
                altitude==0?spp::app::Color{0.31f,0.34f,0.39f,1}:spp::app::Color{0.13f,0.16f,0.20f,1});
        }
        for(int hour=0;hour<=24;hour+=3) {
            const float x=day_chart_.x+day_chart_.w*hour/24;
            view_.line(x,day_chart_.y,x,bottom,{0.10f,0.13f,0.17f,1});
            view_.text(x-9,bottom+3,std::format("{:02d}",hour),0.57f,muted);
        }
        for(std::size_t i=1;i<daily_->altitude.size();++i) {
            const auto n=static_cast<float>(daily_->altitude.size()-1);
            view_.line(day_chart_.x+day_chart_.w*static_cast<float>(i-1)/n,level(daily_->altitude[i-1]),
                day_chart_.x+day_chart_.w*static_cast<float>(i)/n,level(daily_->altitude[i]),accent,1.8f);
        }
        for(const auto& e:h.events) {
            const float x=day_chart_.x+day_chart_.w*static_cast<float>(days(e.time,daily_->anchor));
            const bool crossing=e.kind==spp::moon::HorizonKind::Rise || e.kind==spp::moon::HorizonKind::Set;
            view_.line(x,day_chart_.y-5,x,bottom+2,crossing?ink:muted);
            if(Rect{x-9,day_chart_.y-6,18,day_chart_.h+9}.contains(mx_,my_)) {
                view_.quad({24,top+25,width_-48,24},{0.035f,0.043f,0.057f,1});
                view_.text(24,top+29,std::format("{} {} UTC | center {:+.3f} deg | az {:.3f} deg",
                    spp::moon::name(e.kind),stamp(e.time),e.center_altitude_deg,e.azimuth_deg),0.62f);
            }
            const char* label=e.kind==spp::moon::HorizonKind::Rise?"Rise":e.kind==spp::moon::HorizonKind::Set?"Set":
                e.kind==spp::moon::HorizonKind::UpperTransit?"Upper":e.kind==spp::moon::HorizonKind::LowerTransit?"Lower":"Graze";
            if(button({std::clamp(x-34,24.0f,width_-95),bottom+26,70,24},label)) {
                config.time=e.time;dirty_=true;playing=false;
            }
        }
        const float current=day_chart_.x+day_chart_.w*static_cast<float>(days(config.time,daily_->anchor));
        view_.line(current,day_chart_.y-9,current,bottom+3,ink,1.6f);
    }
    void calculate() {
        if(!dirty_)return;
        const auto before=Clock::now();
        try {
            const spp::Atmosphere atm{.refraction=config.refraction};
            if(config.numerical)observation_=spp::moon::observe(config.time,config.observer,*config.numerical,atm);
            else observation_=config.direct?spp::moon::observe(config.time,config.observer,direct_,atm)
                :spp::moon::observe(config.time,config.observer,cached_,atm);
            valid_=true;error_.clear();
        } catch(const std::exception& e){error_=e.what();valid_=false;playing=false;}
        calculation_us_=std::chrono::duration<double,std::micro>(Clock::now()-before).count();dirty_=false;
    }
    void accept() {
        try {
            if(edit_==0)config.time=timestamp(buffer_);
            else if(edit_==1) {
                const double x=number(buffer_);if(std::abs(x)>90)throw std::invalid_argument("Latitude must be -90 to 90 degrees.");
                config.observer.lat_deg=x;
            } else if(edit_==2)config.observer.lon_deg=number(buffer_);
            else if(edit_==3) {
                const double x=number(buffer_);if(x< -12000 || x>1e7)throw std::invalid_argument("Elevation must be -12000 to 10000000 metres.");
                config.observer.elev_m=x;
            }
            edit_=-1;dirty_=true;error_.clear();
        } catch(const std::exception& e){error_=e.what();}
    }
public:
    Options config;bool playing=false;double speed=1; // simulated days per real second
    Application(Options o,bool background):cached_(o.fit),direct_(o.fit),config(std::move(o)) {
        if(background) {
            worker_=std::make_unique<Worker<Calendar>>([](Calendar& out,spp::app::WorkerCancellation cancel){fill_calendar(out,cancel);},glfwPostEmptyEvent);
            daily_worker_=std::make_unique<Worker<DailySky>>([](DailySky& out,spp::app::WorkerCancellation cancel){fill_calendar(out,cancel);},glfwPostEmptyEvent);
        }
    }
    void advance(double d){config.time=config.time.shifted(d);dirty_=true;}
    void mouse(double x,double y,bool press,bool down){mx_=x;my_=y;pressed_=pressed_||press;down_=down;}
    void character(unsigned int c){if(edit_>=0 && c>=32 && c<127){if(select_all_){buffer_.clear();select_all_=false;}if(buffer_.size()<64)buffer_+=static_cast<char>(c);}}
    void key(GLFWwindow* window,int key,int mods) {
        if(edit_>=0) {
            if(key==GLFW_KEY_ESCAPE){edit_=-1;error_.clear();}
            else if(key==GLFW_KEY_ENTER || key==GLFW_KEY_KP_ENTER)accept();
            else if(key==GLFW_KEY_BACKSPACE){if(select_all_)buffer_.clear();else if(!buffer_.empty())buffer_.pop_back();select_all_=false;}
            else if(key==GLFW_KEY_A && (mods&GLFW_MOD_CONTROL))select_all_=true;
            else if(key==GLFW_KEY_V && (mods&GLFW_MOD_CONTROL)) {
                if(const auto* text=glfwGetClipboardString(window))for(const unsigned char* p=reinterpret_cast<const unsigned char*>(text);*p;++p)character(*p);
            }
            return;
        }
        switch(key) {
        case GLFW_KEY_ESCAPE:glfwSetWindowShouldClose(window,GLFW_TRUE);break;
        case GLFW_KEY_SPACE:playing=!playing;break;
        case GLFW_KEY_LEFT:advance((mods&GLFW_MOD_SHIFT)?-1:-1.0/24);break;
        case GLFW_KEY_RIGHT:advance((mods&GLFW_MOD_SHIFT)?1:1.0/24);break;
        case GLFW_KEY_EQUAL:case GLFW_KEY_KP_ADD:speed=std::min(32.0,speed*2);break;
        case GLFW_KEY_MINUS:case GLFW_KEY_KP_SUBTRACT:speed=std::max(1.0/64,speed/2);break;
        case GLFW_KEY_N:config.time=now();dirty_=true;break;
        case GLFW_KEY_O:north_up_=!north_up_;break;
        case GLFW_KEY_E:earthshine_=!earthshine_;break;
        case GLFW_KEY_G:area_=!area_;break;
        case GLFW_KEY_H:config.events=!config.events;break;
        case GLFW_KEY_L:config.grid=!config.grid;break;
        case GLFW_KEY_F:change_fit();break;
        case GLFW_KEY_D:if(!config.numerical){config.direct=!config.direct;dirty_=true;}break;
        case GLFW_KEY_R:config.refraction=!config.refraction;dirty_=true;break;
        default:break;
        }
    }
    void prime(){calculate();update_calendar();update_daily();}
    double checksum() const {
        if(!valid_)throw std::runtime_error("invalid observation in rendering benchmark: "+error_);
        double sum=observation_.phase.illuminated_fraction+observation_.horizontal.alt_deg+observation_.horizontal.az_deg;
        const auto& surface=observation_.surface;
        sum+=surface.subobserver_lon_deg+surface.subobserver_lat_deg+surface.subsolar_lon_deg+surface.subsolar_lat_deg+surface.pole_position_angle_deg;
        for(const auto& row:surface.body_to_view)for(double x:row)sum+=x;
        if(config.events) {
            if(!daily_ || !daily_->error.empty())throw std::runtime_error("invalid daily sky in rendering benchmark: "+(daily_?daily_->error:"missing result"));
            sum+=daily_->horizon.duration_hours;
            for(double altitude:daily_->altitude)sum+=altitude;
            for(const auto& e:daily_->horizon.events)
                sum+=days(e.time,daily_->anchor)+e.center_altitude_deg+e.azimuth_deg+static_cast<int>(e.kind);
        } else {
            if(!calendar_ || !calendar_->error.empty())throw std::runtime_error("invalid phase calendar in rendering benchmark: "+(calendar_?calendar_->error:"missing result"));
            for(double fraction:calendar_->fraction)sum+=fraction;
            for(const auto& e:calendar_->events)sum+=days(e.time,calendar_->anchor)+static_cast<int>(e.phase);
        }
        return sum;
    }
    std::uint64_t evaluations() const{return config.numerical?0:config.direct?direct_.calls:cached_.source_evaluations();}
    auto schedule_updates() const{return std::pair{calendar_updates_,daily_updates_};}
    void frame(int w,int h) {
        width_=static_cast<float>(w);height_=static_cast<float>(h);view_.begin(w,h);
        if(!config.events && (down_ || pressed_) && requested_anchor_ && timeline_.contains(mx_,my_) && edit_<0) {
            config.time=requested_anchor_->shifted(std::clamp((mx_-static_cast<double>(timeline_.x))
                /static_cast<double>(timeline_.w),0.0,1.0)*42);
            dirty_=true;playing=false;pressed_=false;
        }
        if(config.events && (down_ || pressed_) && requested_day_ && day_chart_.contains(mx_,my_) && edit_<0) {
            config.time=requested_day_->anchor.shifted(std::clamp((mx_-static_cast<double>(day_chart_.x))
                /static_cast<double>(day_chart_.w),0.0,std::nextafter(1.0,0.0)));
            dirty_=true;playing=false;pressed_=false;
        }
        view_.text(24,15,"Moon",1.30f);
        field(0,{155,19,265,32},stamp(config.time));view_.text(428,22,"UTC",0.75f,muted);
        if(button({477,19,68,32},"Now")){config.time=now();dirty_=true;}
        if(button({width_-385,19,160,32},config.numerical?"Numerical Moon":config.fit==Fit::DE405?"DE405/DE406":"LLR"))change_fit();
        if(button({width_-214,19,190,32},config.numerical?"Stored trajectory":config.direct?"Direct coordinates":"Cached coordinates",!config.direct)
            && !config.numerical){config.direct=!config.direct;dirty_=true;}
        if(button({24,73,60,32},"-1d"))advance(-1);
        if(button({91,73,60,32},"-1h"))advance(-1.0/24);
        if(button({158,73,83,32},playing?"Pause":"Play",playing))playing=!playing;
        if(button({248,73,60,32},"+1h"))advance(1.0/24);
        if(button({315,73,60,32},"+1d"))advance(1);
        if(button({390,73,38,32},"-"))speed=std::max(1.0/64,speed/2);
        view_.text(437,78,std::format("{:g} d/s",speed),0.75f,muted);
        if(button({544,73,38,32},"+"))speed=std::min(32.0,speed*2);
        const float side=width_-383;
        view_.text(side,57,"LAT / LON / ELEVATION (m)",0.62f,muted);
        field(1,{side,79,105,32},std::format("{:.4f}",config.observer.lat_deg));
        field(2,{side+112,79,118,32},std::format("{:.4f}",config.observer.lon_deg));
        field(3,{side+237,79,122,32},std::format("{:.1f}",config.observer.elev_m));
        view_.line(24,126,width_-24,126,{0.16f,0.19f,0.23f,1});
        if(button({24,145,158,30},north_up_?"North up":"Zenith up"))north_up_=!north_up_;
        if(button({190,145,151,30},area_?"Lit area":"Surface light",area_))area_=!area_;
        if(button({349,145,132,30},"Earthshine",earthshine_))earthshine_=!earthshine_;
        if(button({489,145,173,30},"Libration grid",config.grid))config.grid=!config.grid;
        calculate();update_calendar();update_daily();
        if(valid_) {
            const auto& o=observation_;const auto& disk=o.disk;
            const float radius=std::max(30.0f,std::min((side-100)*0.5f,(height_-480)*0.5f));
            const float cx=(side+24)*0.5f,cy=200+radius;
            const bool north=north_up_ || !disk.zenith_defined;
            const double chi=(disk.bright_limb_deg-(north?0:disk.parallactic_angle_deg))*spp::deg2rad;
            const double phase=disk.phase_angle_deg*spp::deg2rad;
            view_.disk(cx,cy,radius,{static_cast<float>(-std::sin(phase)*std::sin(chi)),
                static_cast<float>(std::sin(phase)*std::cos(chi)),static_cast<float>(std::cos(phase))},earthshine_,area_);
            if(config.grid)draw_grid(view_,cx,cy,radius,o.surface,north?0:disk.parallactic_angle_deg*spp::deg2rad);
            view_.text(cx-18,cy-radius-27,north?"N":"UP",0.7f,muted);
            view_.text(cx-radius-27,cy-13,north?"E":"L",0.7f,muted);
            const char* orientation=!north_up_ && !disk.zenith_defined?"Zenith direction undefined; north shown":
                north?"Celestial north up, east left":"Zenith up, before atmospheric distortion";
            view_.text(cx-175,cy+radius+12,config.grid?(north?"IAU2009 ME grid | Celestial north up":"IAU2009 ME grid | Zenith up"):orientation,0.68f,muted);
            view_.text(side,153,o.phase.name(),1.06f);
            view_.text(side,190,std::format("{:.2f}%",100*disk.illuminated_fraction),1.95f);
            view_.text(side,245,"of the visible disk illuminated",0.70f,muted);
            float y=294;
            auto row=[&](std::string_view label,const std::string& value) {view_.text(side,y,label,0.76f,muted);view_.text(side+164,y,value,0.76f);y+=26;};
            row("Altitude",std::format("{:+.4f} deg",o.horizontal.alt_deg));
            row("Azimuth",std::format("{:.4f} deg",o.horizontal.az_deg));
            row("Distance",std::format("{:.0f} km",o.distance_km));
            row("Diameter",std::format("{:.3f} arcmin",o.semidiameter_deg*120));
            row("Phase angle",std::format("{:.4f} deg",disk.phase_angle_deg));
            row("Elongation",std::format("{:.4f} deg",o.phase.elongation_deg));
            if(config.grid) {
                row("Moon pole PA",o.surface.pole_defined?std::format("{:.3f} deg",o.surface.pole_position_angle_deg):"undefined");
                row("Sun ME lon",std::format("{:+.4f} deg",o.surface.subsolar_lon_deg));
                row("ME lon (east)",std::format("{:+.4f} deg",o.surface.subobserver_lon_deg));
                row("ME latitude",std::format("{:+.4f} deg",o.surface.subobserver_lat_deg));
            } else {
            row("Bright limb",disk.bright_limb_defined?std::format("{:.3f} deg",disk.bright_limb_deg):"undefined");
            row("Parallactic",disk.zenith_defined?std::format("{:+.3f} deg",disk.parallactic_angle_deg):"undefined");
            if(calendar_ && calendar_->error.empty()) {
                std::optional<Time> previous,next;
                for(const auto& e:calendar_->events)if(e.phase==spp::moon::PrincipalPhase::NewMoon) {
                    if(days(e.time,config.time)<=0)previous=e.time;else if(!next)next=e.time;
                }
                if(previous)row("Age",std::format("{:.4f} days",days(config.time,*previous)));
                if(previous && next)row("Lunation",std::format("{:.4f} days",days(*next,*previous)));
            }
            }
            if(button({side,y+12,359,30},config.refraction?"Atmospheric refraction on":"Atmospheric refraction off",config.refraction)){config.refraction=!config.refraction;dirty_=true;}
        }
        if(!error_.empty())view_.text(24,190,error_,0.76f,{1,0.50f,0.40f,1});
        const float top=height_-211;
        view_.line(24,top-12,width_-24,top-12,{0.16f,0.19f,0.23f,1});
        if(button({24,top-4,145,28},"Phases",!config.events))config.events=false;
        if(button({177,top-4,187,28},"Today's sky",config.events)){config.events=true;update_daily();}
        if(config.events)daily_panel(top);
        else {
        timeline_={24,top+34,width_-48,118};
        const float bottom=timeline_.y+77;
        for(int i=0;i<=6;++i) {
            const float x=timeline_.x+timeline_.w*static_cast<float>(i)/6;
            view_.line(x,timeline_.y,x,bottom,{0.10f,0.13f,0.17f,1});
            if(requested_anchor_)view_.text(std::clamp(x-40,24.0f,width_-125),bottom+33,spp::iso_date(requested_anchor_->shifted(i*7).to_civil().date()),0.58f,muted);
        }
        if(calendar_ && calendar_->error.empty() && !calendar_->fraction.empty()) {
            for(std::size_t i=1;i<calendar_->fraction.size();++i) {
                const float a=static_cast<float>(i-1)/static_cast<float>(calendar_->fraction.size()-1),b=static_cast<float>(i)/static_cast<float>(calendar_->fraction.size()-1);
                view_.line(timeline_.x+timeline_.w*a,bottom-static_cast<float>(calendar_->fraction[i-1])*72,
                    timeline_.x+timeline_.w*b,bottom-static_cast<float>(calendar_->fraction[i])*72,accent,1.8f);
            }
            for(const auto& e:calendar_->events) {
                const double d=days(e.time,calendar_->anchor);if(d<0 || d>=42)continue;
                const float x=timeline_.x+timeline_.w*static_cast<float>(d/42);
                view_.line(x,timeline_.y-5,x,bottom+3,{0.36f,0.39f,0.43f,1});
                const char* label=e.phase==spp::moon::PrincipalPhase::NewMoon?"New":e.phase==spp::moon::PrincipalPhase::FullMoon?"Full":
                    e.phase==spp::moon::PrincipalPhase::FirstQuarter?"First":"Last";
                if(Rect{x-8,timeline_.y-5,16,85}.contains(mx_,my_))
                    view_.text(std::clamp(x-130,24.0f,width_-275),top-30,stamp(e.time)+" UTC",0.72f,ink);
                if(button({std::clamp(x-28,24.0f,width_-82),bottom+5,58,24},label)){config.time=e.time;dirty_=true;playing=false;}
            }
        } else view_.text(24,timeline_.y+30,calendar_?calendar_->error:"Computing lunar calendar...",0.78f,muted);
        if(requested_anchor_) {
            const float x=timeline_.x+timeline_.w*static_cast<float>(days(config.time,*requested_anchor_)/42);
            view_.line(x,timeline_.y-9,x,bottom+4,ink,1.5f);
        }
        }
        view_.text(24,height_-39,"Space: play   Arrows: hour   H: daily sky   L: libration grid   Drag chart: time",0.62f,muted);
        view_.text(width_-370,height_-39,std::format("Query {:8.1f} us   {}",calculation_us_,config.numerical?"integrated":config.direct?"direct":"cached"),0.70f,muted);
        view_.end();pressed_=false;
    }
};

void capture(const std::string& path,int w,int h) {
    std::vector<unsigned char> pixels(static_cast<std::size_t>(w)*static_cast<std::size_t>(h)*3);
    glPixelStorei(GL_PACK_ALIGNMENT,1);glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,pixels.data());
    std::ofstream file(path,std::ios::binary);if(!file)throw std::runtime_error("cannot create capture: "+path);
    file<<"P6\n"<<w<<' '<<h<<"\n255\n";
    for(int row=h-1;row>=0;--row)file.write(reinterpret_cast<const char*>(pixels.data()+static_cast<std::size_t>(row)*static_cast<std::size_t>(w)*3),w*3);
    if(!file)throw std::runtime_error("cannot write capture: "+path);
}
void renderer_info() {
    std::println("OpenGL {} | {} | {}",reinterpret_cast<const char*>(glGetString(GL_VERSION)),
        reinterpret_cast<const char*>(glGetString(GL_VENDOR)),reinterpret_cast<const char*>(glGetString(GL_RENDERER)));
}
void initialize_gl() {
    glewExperimental=GL_TRUE;const GLenum status=glewInit();
    if(status!=GLEW_OK && status!=GLEW_ERROR_NO_GLX_DISPLAY)
        throw std::runtime_error(std::string("GLEW: ")+reinterpret_cast<const char*>(glewGetErrorString(status)));
    while(glGetError()!=GL_NO_ERROR){} // GLEW probes an extension string on core contexts.
    if(!glCreateShader || !glGenVertexArrays)throw std::runtime_error("OpenGL 3.3 entry points unavailable");
    renderer_info();
}
int benchmark(const Options& o) {
    Application app(o,false);app.prime();app.frame(1280,900);glFinish();
    std::vector<double> times;times.reserve(static_cast<std::size_t>(o.frames));double checksum=0;
    const auto evaluations=app.evaluations();
    const auto updates=app.schedule_updates();
    for(int i=0;i<o.frames;++i) {
        app.advance(o.step_seconds/spp::seconds_per_day);
        const auto before=Clock::now();app.frame(1280,900);glFinish();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-before).count());checksum+=app.checksum();
        if(glGetError()!=GL_NO_ERROR)throw std::runtime_error("OpenGL error during rendering benchmark");
    }
    if(!std::isfinite(checksum))throw std::runtime_error("non-finite rendering benchmark checksum");
    if(!o.capture.empty())capture(o.capture,1280,900);
    std::sort(times.begin(),times.end());double total=0;for(double t:times)total+=t;
    std::println("Moon render + observation: {} frames at 1280x900, {:.3f} s step/frame, {}",o.frames,o.step_seconds,o.numerical?"integrated":o.direct?"direct":"cached");
    std::println("  mean {:.3f} ms, median {:.3f} ms, p95 {:.3f} ms, max {:.3f} ms; {:.1f} frames/s",
        total/o.frames,times[times.size()/2],times[std::min(times.size()-1,times.size()*95/100)],times.back(),1000*o.frames/total);
    std::println("  {} observation ELP {} evaluations, checksum {:.9f}; initial calendar/warmup excluded, glFinish included",
        app.evaluations()-evaluations,o.direct?"position":"endpoint state",checksum);
    std::println("  Calendar refreshes inside the measured timeline are included; no swap/vsync wait.");
    std::println("  Panel: {}; measured refreshes: {} phase calendars, {} daily skies ({} schedule evaluation)",
        o.events?"daily sky":"phases",app.schedule_updates().first-updates.first,app.schedule_updates().second-updates.second,
        o.numerical?"integrated":"cached");
    return 0;
}
}  // namespace

int main(int argc,char** argv) {
    try {
        const auto o=options(argc,argv);
        if(o.offscreen){spp::app::Offscreen context(1280,900);initialize_gl();return benchmark(o);}
        glfwSetErrorCallback([](int,const char* text){std::fprintf(stderr,"GLFW: %s\n",text);});
        if(!glfwInit())throw std::runtime_error("GLFW needs a desktop display; use --offscreen --benchmark=N for headless performance runs.");
        struct GLFWCleanup{~GLFWCleanup(){glfwTerminate();}} cleanup;
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
        glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);glfwWindowHint(GLFW_SAMPLES,4);
#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT,GLFW_TRUE);
#endif
        if(o.frames)glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
        GLFWwindow* window=glfwCreateWindow(1280,900,"SPP - Moon",nullptr,nullptr);
        if(!window)throw std::runtime_error("cannot create the OpenGL 3.3 Moon window");
        struct WindowCleanup{GLFWwindow* window;~WindowCleanup(){glfwDestroyWindow(window);}} window_cleanup{window};
        glfwMakeContextCurrent(window);initialize_gl();glfwSwapInterval(o.frames?0:1);
        if(o.frames)return benchmark(o);
        glfwSetWindowSizeLimits(window,1100,820,GLFW_DONT_CARE,GLFW_DONT_CARE);
        Application app(o,true);glfwSetWindowUserPointer(window,&app);
        glfwSetKeyCallback(window,[](GLFWwindow* w,int key,int,int action,int mods){
            if(action==GLFW_PRESS || (action==GLFW_REPEAT && (key==GLFW_KEY_LEFT || key==GLFW_KEY_RIGHT
                || key==GLFW_KEY_BACKSPACE || key==GLFW_KEY_EQUAL || key==GLFW_KEY_MINUS)))
                static_cast<Application*>(glfwGetWindowUserPointer(w))->key(w,key,mods);
        });
        glfwSetCharCallback(window,[](GLFWwindow* w,unsigned int c){static_cast<Application*>(glfwGetWindowUserPointer(w))->character(c);});
        glfwSetMouseButtonCallback(window,[](GLFWwindow* w,int button,int action,int) {
            if(button!=GLFW_MOUSE_BUTTON_LEFT)return;
            double x=0,y=0;int fw=0,fh=0,ww=0,wh=0;
            glfwGetCursorPos(w,&x,&y);glfwGetFramebufferSize(w,&fw,&fh);glfwGetWindowSize(w,&ww,&wh);
            if(ww && wh)static_cast<Application*>(glfwGetWindowUserPointer(w))->mouse(x*fw/ww,y*fh/wh,action==GLFW_PRESS,action==GLFW_PRESS);
        });
        bool was_down=false;auto last=Clock::now();
        while(!glfwWindowShouldClose(window)) {
            glfwPollEvents();const auto current=Clock::now();
            const double elapsed=std::min(0.1,std::chrono::duration<double>(current-last).count());last=current;
            if(app.playing)app.advance(elapsed*app.speed);
            int w=0,h=0,ww=0,wh=0;glfwGetFramebufferSize(window,&w,&h);glfwGetWindowSize(window,&ww,&wh);
            if(w==0 || h==0 || ww==0 || wh==0){glfwWaitEventsTimeout(0.25);continue;}
            double x=0,y=0;glfwGetCursorPos(window,&x,&y);
            const bool down=glfwGetMouseButton(window,GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS;
            app.mouse(x*w/ww,y*h/wh,down&&!was_down,down);was_down=down;
            app.frame(w,h);glfwSwapBuffers(window);
            if(!app.playing)glfwWaitEvents();
        }
        return 0;
    } catch(const std::exception& e){std::fprintf(stderr,"spp_moon: %s\n",e.what());return 1;}
}
