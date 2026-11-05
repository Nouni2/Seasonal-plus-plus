// Cache-residency and eviction gates coupled to Earth-cache timings.
#include "earth_state.hpp"

#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using spp::detail::EarthState;
using Epoch = spp::Instant<spp::TT>;

std::array<double,5> values(const EarthState& s) {
    return {s.l_rad,s.b_rad,s.r_au,s.dpsi_arcsec,s.deps_arcsec};
}
void equal(const EarthState& a,const EarthState& b) {
    const auto x=values(a),y=values(b);
    for(std::size_t i=0;i<x.size();++i)
        if(!std::isfinite(x[i]) || std::bit_cast<std::uint64_t>(x[i])!=std::bit_cast<std::uint64_t>(y[i]))
            throw std::runtime_error("cache eviction changed an Earth-state value");
}
void workload(const char* name,std::int64_t stride) {
    spp::detail::clear_earth_cache();
    const auto capacity=spp::detail::earth_cache_stats().capacity;
    std::vector<Epoch> epochs;
    std::vector<EarthState> expected;
    const auto begin=Clock::now();
    for(std::size_t i=0;i<capacity;++i) {
        epochs.push_back({spp::jd_j2000_day+(static_cast<std::int64_t>(i)-128)*stride*16,0.375});
        expected.push_back(spp::detail::earth_state(epochs.back()));
    }
    const auto cold=std::chrono::duration<double,std::milli>(Clock::now()-begin).count();
    const auto before=spp::detail::earth_cache_stats();
    const auto start=Clock::now();
    double checksum=0;
    for(int round=0;round<20;++round)for(std::size_t i=0;i<epochs.size();++i) {
        const auto s=spp::detail::earth_state(epochs[i]);
        equal(s,expected[i]);
        checksum+=s.l_rad+s.b_rad+s.r_au+s.dpsi_arcsec+s.deps_arcsec;
    }
    const auto us=std::chrono::duration<double,std::micro>(Clock::now()-start).count()/(20*epochs.size());
    const auto after=spp::detail::earth_cache_stats();
    if(before.builds!=capacity || before.builds!=after.builds)
        throw std::runtime_error("resident Earth segments were unnecessarily rebuilt");
    // Repeated eviction and hash-index rebuilds must preserve every value.
    for(int round=0;round<4;++round) {
        for(std::size_t i=0;i<capacity*2;++i)
            (void)spp::detail::earth_state({spp::jd_j2000_day+(static_cast<std::int64_t>(i)+4096+round*1024)*16,0.8125});
        for(std::size_t i=epochs.size();i-- >0;)equal(spp::detail::earth_state(epochs[i]),expected[i]);
    }
    // An equivalent noncanonical split must address the same segment/state.
    for(std::size_t i=0;i<epochs.size();++i)
        equal(spp::detail::earth_state({epochs[i].day-1,epochs[i].frac+1}),expected[i]);
    std::printf("  %-24s %zu segments: cold %.3f ms, resident %.3f us/query, %llu rebuilds, checksum %.12g\n",
        name,capacity,cold,us,static_cast<unsigned long long>(after.builds-before.builds),checksum);
}
}

int main() {
    try {
        const auto s=spp::detail::earth_cache_stats();
        std::printf("Earth cache performance: %zu slots, %zu bytes/thread; resident time includes bitwise value gates\n",s.capacity,s.bytes);
        workload("consecutive segments",1);
        workload("old modulo collisions",64);
        std::puts("PASS: Earth residency, eviction, split-time and numerical identity gates");
        return 0;
    } catch(const std::exception& e) {std::fprintf(stderr,"FAIL Earth cache benchmark: %s\n",e.what());return 1;}
}
