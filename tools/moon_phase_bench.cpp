// Matched-time independent phase gates coupled to direct/cached search cost.
#include "spp/spp.hpp"
#include "spp/model/elp_mpp02.hpp"
#include "error_stats.hpp"
#include "phase_reference.hpp"
#include "core/aberration.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <format>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Clock=std::chrono::steady_clock;
using spp::check::phase_reference::Arc;
using spp::check::phase_reference::read;
using spp::check::phase_reference::event;
using spp::check::phase_reference::seconds;
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}

void aberration_benchmark() {
    using namespace spp::detail;
    struct Case {Vector n,velocity;};
    std::vector<Case> cases;
    double inverse_error=0,length_error=0;
    require(norm(aberration({0,1,0},light_speed_km_day*Vector{0.6,0,0})-Vector{0.6,0.8,0})<1e-14,
        "relativistic transverse-ray reference failed");
    for(int j=0;j<24;++j)for(double latitude:{-1.4,-0.4,0.4,1.4})for(double beta:{0.0,1e-4,0.01,0.6}) {
        const double l=j*(2*spp::pi)/24;
        const Vector n{std::cos(latitude)*std::cos(l),std::cos(latitude)*std::sin(l),std::sin(latitude)};
        const Vector v=light_speed_km_day*Vector{0.6*beta,-0.8*beta,0};
        const auto apparent=aberration(n,v);
        inverse_error=std::max(inverse_error,norm(aberration(apparent,(-1.0)*v)-n));
        length_error=std::max(length_error,std::abs(norm(apparent)-1));
        cases.push_back({n,v});
    }
    double checksum=0;
    const auto before=Clock::now();
    for(int repeat=0;repeat<256;++repeat)for(const auto& c:cases) {
        const auto v=aberration(c.n,c.velocity);
        checksum+=v[0]+v[1]+v[2];
    }
    const double us=std::chrono::duration<double,std::micro>(Clock::now()-before).count()/(256*cases.size());
    std::printf("Illumination aberration: %.4f us/query; inverse error %.3g, length error %.3g; checksum %.9g\n",
        us,inverse_error,length_error,checksum);
    require(inverse_error<1e-14 && length_error<1e-14 && std::isfinite(checksum),"illumination aberration numerical gate failed");
}

void benchmark(const std::vector<Arc>& arcs) {
    spp::model::elp_mpp02::Ephemeris direct;
    spp::model::elp_mpp02::CachedEphemeris cached;
    spp::check::ErrorStats physical,interpolation;
    std::vector<spp::Time> exact,fast;
    const auto start=Clock::now();
    for(const auto& arc:arcs)exact.push_back(event(arc,direct));
    const double direct_ms=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
    const auto cold=Clock::now();
    for(const auto& arc:arcs)fast.push_back(event(arc,cached));
    const double cold_ms=std::chrono::duration<double,std::milli>(Clock::now()-cold).count();
    for(std::size_t i=0;i<arcs.size();++i) {
        const auto label=std::format("{} {}",spp::moon::name(static_cast<spp::moon::PrincipalPhase>(arcs[i].kind)),exact[i].iso_utc());
        physical.add(seconds(exact[i].tt(),arcs[i].root),label);
        interpolation.add(seconds(fast[i].tt(),exact[i].tt()),label);
    }
    const auto calls=cached.source_evaluations();
    std::vector<double> times;
    double checksum=0;
    for(int round=0;round<5;++round) {
        const auto before=Clock::now();
        for(const auto& arc:arcs)checksum+=seconds(event(arc,cached).tt(),arc.root);
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-before).count());
    }
    std::sort(times.begin(),times.end());
    std::printf("Phase short-arc searches: direct %.3f ms; initially empty lunar cache %.3f ms; resident median %.3f ms / %zu events\n",
        direct_ms,cold_ms,times[times.size()/2],arcs.size());
    std::printf("  %llu cached endpoint states to fill; %llu in five repeat batches; checksum %.9f s\n",
        static_cast<unsigned long long>(calls),static_cast<unsigned long long>(cached.source_evaluations()-calls),checksum);
    const bool independent=physical.report("phase vs DE441 (TT)","s",0.5);
    const bool adapter=interpolation.report("cached vs direct time","s",0.0002);
    require(independent && adapter,"matched-time phase accuracy gate failed");
    require(std::isfinite(checksum) && cached.source_evaluations()==calls,"resident phase workload unnecessarily refilled its cache");
}
}

int main(int argc,char** argv) {
    try {
        require(argc==2,"usage: spp_moon_phase_bench moon_horizons_phase_tt.csv");
        aberration_benchmark();benchmark(read(argv[1]));
        std::puts("PASS: matched-TT independent phase gates and performance");return 0;
    } catch(const std::exception& e){std::fprintf(stderr,"FAIL phase benchmark: %s\n",e.what());return 1;}
}
