#pragma once

#include "spp/spp.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Shared independent short-arc reference reader. It loads no Spp ephemeris;
// TT roots come only from the pinned Horizons samples and an unused control.
namespace spp::check::phase_reference {
using TT=spp::Instant<spp::TT>;
inline void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
inline double seconds(TT a,TT b){return (static_cast<double>(a.day-b.day)+(a.frac-b.frac))*86400;}
inline TT split(double jd){const auto day=static_cast<std::int64_t>(std::floor(jd));return {day,jd-static_cast<double>(day)};}

inline spp::Time at_tt(TT target) {
    auto t=spp::Time::from_jd(target.day,target.frac);
    for(int i=0;i<4;++i)t=t.shifted(-seconds(t.tt(),target)/86400);
    require(std::abs(seconds(t.tt(),target))<1e-6,"TT-to-UTC inversion failed");
    return t;
}
struct Arc {
    int kind=0;
    std::array<double,4> jd{},longitude{};
    spp::Time start,end;
    TT root;
};

inline std::vector<Arc> read(const char* path) {
    std::ifstream file(path);require(file.good(),"cannot open JPL TT phase fixture");
    struct Samples {int kind;std::vector<double> jd,longitude;};
    std::vector<Samples> input;
    for(std::string line;std::getline(file,line);) {
        if(line.empty() || line.front()=='#')continue;
        std::vector<std::string> c;std::istringstream row(line);
        for(std::string x;std::getline(row,x,',');)c.push_back(x);
        require(c.size()==5,"invalid JPL phase fixture row");
        const auto i=static_cast<std::size_t>(std::stoul(c[0]));const int kind=std::stoi(c[1]);
        require(kind>=0 && kind<=3 && i<=input.size(),"invalid phase fixture group");
        if(i==input.size())input.push_back({kind,{},{}});
        require(input[i].kind==kind,"phase kind changed within one short arc");
        input[i].jd.push_back(std::stod(c[2]));
        input[i].longitude.push_back(spp::wrap180(std::stod(c[3])-std::stod(c[4])-90*kind));
    }
    std::vector<Arc> result;
    double control_error=0;
    for(const auto& c:input) {
        require(c.jd.size()==4,"phase arc needs four independent samples");
        std::array<double,4> x{};
        for(std::size_t i=0;i<4;++i) {
            x[i]=(c.jd[i]-c.jd[1])*86400;
            require(std::isfinite(c.longitude[i]) && (i==0 || c.jd[i]>c.jd[i-1]),"nonfinite or unordered phase reference");
        }
        // Fit three reference samples and check against the fourth, which
        // is not used in the root. Quantized Horizons longitudes are 1e-7 deg.
        auto reference=[&](double t) {
            double y=0;
            for(int i:{0,1,3}) {
                double term=c.longitude[static_cast<std::size_t>(i)];
                for(int j:{0,1,3})if(i!=j)term*=(t-x[static_cast<std::size_t>(j)])/(x[static_cast<std::size_t>(i)]-x[static_cast<std::size_t>(j)]);
                y+=term;
            }
            return y;
        };
        control_error=std::max(control_error,std::abs(reference(x[2])-c.longitude[2]));
        require(c.longitude[0]<0 && c.longitude[3]>0,"JPL phase root not bracketed");
        double a=x[0],b=x[3];
        while(b-a>1e-5) {const double m=(a+b)/2;if(reference(m)<0)a=m;else b=m;}
        auto root=split(c.jd[1]);root.frac+=(a+b)/(2*86400);
        if(root.frac<0){--root.day;root.frac+=1;}else if(root.frac>=1){++root.day;root.frac-=1;}
        Arc arc{c.kind,{},{},at_tt(split(c.jd[0])),at_tt(split(c.jd[3])),root};
        std::copy(c.jd.begin(),c.jd.end(),arc.jd.begin());
        std::copy(c.longitude.begin(),c.longitude.end(),arc.longitude.begin());
        result.push_back(arc);
    }
    require(result.size()==62,"JPL phase coverage changed; review the fixture");
    std::printf("Independent DE441 TT phase arcs: %zu events; unused-sample quadratic residual %.9g arcsec\n",result.size(),control_error*3600);
    require(control_error<3e-7,"JPL short-arc interpolation control failed");
    return result;
}

template<class Model>
spp::Time event(const Arc& arc,Model& model) {
    const auto roots=spp::moon::phases(arc.start,arc.end,model,0.0001);
    require(roots.size()==1 && static_cast<int>(roots.front().phase)==arc.kind,"phase search count/kind differs from JPL");
    return roots.front().time;
}

} // namespace spp::check::phase_reference
