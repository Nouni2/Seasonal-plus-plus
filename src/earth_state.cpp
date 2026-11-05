#include "earth_state.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

#include "core/split.hpp"
#include "spp/core/constants.hpp"
#include "spp/model/iau1980.hpp"
#include "spp/model/vsop87d.hpp"

namespace spp::detail {

EarthState earth_state_direct(double T) {
    const auto e = model::vsop87d::earth(T / 10.0);
    const auto n = model::iau1980::nutation(T);
    return {e.l_rad, e.b_rad, e.r_au, n.dpsi_arcsec, n.deps_arcsec};
}

namespace {

// Segments are kSegmentDays long, aligned on J2000.0 (a TT noon), and each
// holds a degree-kDegree Chebyshev series per quantity. The series fits the
// residual from an offset (the value at the middle node): VSOP longitude
// reaches ~2e4 rad, and fitting it directly puts that magnitude's rounding into
// every coefficient. tools/cache_check.cpp measures the effect: without the
// offset the worst error in L grows 2.6x (1.4 bits).
constexpr std::int64_t kSegmentDays = 16;
constexpr int kDegree = 24;
constexpr int kNodes = kDegree + 1;
constexpr int kQuantities = 5;
constexpr std::size_t kCacheSlots = 256, kBuckets = 512;
constexpr std::uint16_t kTombstone = 0xffff;

struct Segment {
    std::int64_t index = std::numeric_limits<std::int64_t>::min();
    std::array<double, kQuantities> offset{};
    std::array<std::array<double, kQuantities>, kNodes> c{};   // [degree][quantity]
    std::uint64_t stamp = 0;
    std::size_t bucket = 0;
};

thread_local std::array<Segment, kCacheSlots> cache{};
thread_local std::array<std::uint16_t, kBuckets> lookup{};
thread_local std::size_t tombstones = 0;
thread_local std::uint64_t clock = 0, hits = 0, builds = 0;

std::array<double, kQuantities> as_array(const EarthState& s) {
    return {s.l_rad, s.b_rad, s.r_au, s.dpsi_arcsec, s.deps_arcsec};
}

void build(Segment& seg, std::int64_t index) {
    // Sample at the Chebyshev nodes x_k = cos(π (k + 1/2) / n), then project:
    // c_j = (2/n) Σ_k f(x_k) cos(j π (k + 1/2) / n).
    std::array<std::array<double, kNodes>, kQuantities> f{};
    const double start_days = static_cast<double>(index * kSegmentDays);
    const double half = 0.5 * static_cast<double>(kSegmentDays);

    for (int k = 0; k < kNodes; ++k) {
        const double x = std::cos(pi * (k + 0.5) / kNodes);
        const double days = start_days + half * (x + 1.0);
        const auto v = as_array(earth_state_direct(days / days_per_century));
        for (int q = 0; q < kQuantities; ++q) f[q][k] = v[q];
    }

    for (int q = 0; q < kQuantities; ++q) {
        const double offset = f[q][kNodes / 2];
        seg.offset[q] = offset;
        for (int j = 0; j < kNodes; ++j) {
            double sum = 0.0;
            for (int k = 0; k < kNodes; ++k)
                sum += (f[q][k] - offset) * std::cos(pi * j * (k + 0.5) / kNodes);
            seg.c[j][q] = (2.0 / kNodes) * sum;
        }
        seg.c[0][q] *= 0.5;
    }
    seg.index = index;
}

std::int64_t floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    return (a % b != 0 && (a < 0) != (b < 0)) ? q - 1 : q;
}

const Segment& segment(std::int64_t index) {
    auto bucket = [](std::int64_t key) {
        auto hash = static_cast<std::uint64_t>(key);
        hash ^= hash >> 30; hash *= 0xbf58476d1ce4e5b9ULL;
        hash ^= hash >> 27; hash *= 0x94d049bb133111ebULL;
        hash ^= hash >> 31;
        return static_cast<std::size_t>(hash & (kBuckets - 1));
    };
    if (clock == std::numeric_limits<std::uint64_t>::max()) {
        // Preserve recency order across the counter wrap, without changing
        // coefficients or invalidating already resident segments.
        std::array<Segment*, kCacheSlots> ordered{};
        for (std::size_t i = 0; i < cache.size(); ++i) ordered[i] = &cache[i];
        std::sort(ordered.begin(), ordered.end(), [](auto a, auto b){return a->stamp < b->stamp;});
        clock = 0;
        for (auto* s : ordered) if (s->stamp) s->stamp = ++clock;
    }
    ++clock;
    if (tombstones > kBuckets / 4) {
        lookup.fill(0); tombstones = 0;
        for (std::size_t i = 0; i < cache.size(); ++i) if (cache[i].stamp) {
            auto b = bucket(cache[i].index);
            while (lookup[b]) b = (b + 1) & (kBuckets - 1);
            lookup[b] = static_cast<std::uint16_t>(i + 1); cache[i].bucket = b;
        }
    }
    auto pos = bucket(index), insertion = kBuckets;
    for (std::size_t probe = 0; probe < kBuckets; ++probe, pos = (pos + 1) & (kBuckets - 1)) {
        const auto slot = lookup[pos];
        if (!slot) { if (insertion == kBuckets) insertion = pos; break; }
        if (slot == kTombstone) { if (insertion == kBuckets) insertion = pos; continue; }
        auto& s = cache[slot - 1];
        if (s.index == index) { s.stamp = clock; ++hits; return s; }
    }
    if (insertion == kBuckets) throw std::logic_error("Earth cache index invariant failed");
    auto* victim = &cache[0];
    for (auto& s : cache) if (s.stamp < victim->stamp) victim = &s;
    if (victim->stamp) { lookup[victim->bucket] = kTombstone; ++tombstones; }
    build(*victim, index); ++builds;
    if (lookup[insertion] == kTombstone) --tombstones;
    lookup[insertion] = static_cast<std::uint16_t>(victim - cache.data() + 1);
    victim->stamp = clock; victim->bucket = insertion;
    return *victim;
}

EarthState evaluate(Instant<TT> t) {
    // frac is in [0, 1), so the integer day alone selects the segment.
    const std::int64_t days = t.day - jd_j2000_day;
    const std::int64_t index = floor_div(days, kSegmentDays);
    const Segment& seg = segment(index);

    // Local coordinate in [-1, 1), from the split time to keep full precision.
    const double local = static_cast<double>(days - index * kSegmentDays) + t.frac;
    const double x = local * (2.0 / static_cast<double>(kSegmentDays)) - 1.0;

    // Clenshaw recurrence, the five quantities stepped together so their
    // independent dependency chains overlap (and vectorize).
    std::array<double, kQuantities> b1{}, b2{};
    const double two_x = 2.0 * x;
    for (int j = kDegree; j >= 1; --j) {
        for (int q = 0; q < kQuantities; ++q) {
            const double b0 = two_x * b1[q] - b2[q] + seg.c[j][q];
            b2[q] = b1[q];
            b1[q] = b0;
        }
    }
    std::array<double, kQuantities> out{};
    for (int q = 0; q < kQuantities; ++q) out[q] = seg.offset[q] + (x * b1[q] - b2[q] + seg.c[0][q]);
    return {out[0], out[1], out[2], out[3], out[4]};
}

struct Memo {
    std::int64_t day = std::numeric_limits<std::int64_t>::min();
    double frac = 0.0;
    EarthState state{};
};

thread_local Memo last;

}  // namespace

// The pipeline asks for the same instant twice (nutation, then the Sun).
EarthState earth_state(Instant<TT> t) {
    t = normalized(t);
    if (t.day != last.day || t.frac != last.frac) last = {t.day, t.frac, evaluate(t)};
    return last.state;
}

EarthCacheStats earth_cache_stats() {
    return {hits, builds, kCacheSlots, sizeof(cache) + sizeof(lookup) + sizeof(last)
        + sizeof(clock) + sizeof(hits) + sizeof(builds) + sizeof(tombstones)};
}

void clear_earth_cache() {
    for (auto& s : cache) s.stamp = 0;
    lookup.fill(0); tombstones = 0; clock = 0; hits = 0; builds = 0;
    last = Memo{};
}

}  // namespace spp::detail
