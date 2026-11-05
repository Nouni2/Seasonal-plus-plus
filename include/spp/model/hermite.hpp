#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "spp/core/detail/split.hpp"
#include "spp/core/state.hpp"

namespace spp::model {
namespace detail {
template<class State> struct PositionOf;
template<class Origin, class Frame>
struct PositionOf<CartesianState<Origin, Frame>> {
    using type = CartesianPosition<Origin, Frame>;
};
}  // namespace detail

// Bounded cubic Hermite interpolation of an owned geometric state provider.
// Endpoint velocities must be derivatives in km per day of Scale, in the same
// fixed frame as the positions. There are no body-specific equations here.
// A new interval needs two full states; adjacent intervals reuse their shared
// endpoint. A hit evaluates a cubic without calling the source. Accuracy
// depends on the source and interval length and must be checked before using
// this adapter for a different ephemeris.
// Own one instance per thread, or synchronize externally. No allocation/I/O.
template<class Model, class Scale = TDB, int IntervalsPerDay = 288>
class HermiteEphemeris {
    static_assert(IntervalsPerDay > 0 && IntervalsPerDay <= 86400);
public:
    static constexpr int intervals_per_day = IntervalsPerDay;
    using State = std::invoke_result_t<Model&, Instant<Scale>>;
    using Position = typename detail::PositionOf<State>::type;

    HermiteEphemeris(Model source, Instant<Scale> first, Instant<Scale> last)
        : source_(std::move(source)), first_(normalize(first)), last_(normalize(last)) {
        constexpr auto lo = std::numeric_limits<std::int64_t>::min() / IntervalsPerDay + 1;
        constexpr auto hi = std::numeric_limits<std::int64_t>::max() / IntervalsPerDay - 1;
        if (!less(first_, last_) || first_.day < lo || last_.day > hi)
            throw std::invalid_argument("Hermite ephemeris: invalid interpolation domain");
    }

    State operator()(Instant<Scale> t) {
        t = checked(t);
        if (equal(t, last_)) return sample(t);
        const auto& s = segment(t);
        const double u = coordinate(s, t);
        State out{};
        for (std::size_t q = 0; q < 3; ++q) {
            out.position_km[q] = s.c[0][q] + u*(s.c[1][q] + u*(s.c[2][q] + u*s.c[3][q]));
            out.velocity_km_per_day[q] = (s.c[1][q] + u*(2*s.c[2][q] + 3*u*s.c[3][q])) / s.width;
        }
        return out;
    }

    Position position(Instant<Scale> t) {
        t = checked(t);
        if (equal(t, last_)) return {sample(t).position_km};
        const auto& s = segment(t);
        const double u = coordinate(s, t);
        Position out{};
        for (std::size_t q = 0; q < 3; ++q)
            out.position_km[q] = s.c[0][q] + u*(s.c[1][q] + u*(s.c[2][q] + u*s.c[3][q]));
        return out;
    }

    // Prevent changing the source underneath populated segments.
    const Model& source() const { return source_; }
    std::uint64_t source_evaluations() const { return source_evaluations_; }
    void clear() {
        for (auto& s : cache_) s.stamp = 0;
        index_.fill(0);tombstones_=0;
        for (auto& n : endpoints_) n.valid = false;
        clock_ = 0;
        source_evaluations_ = 0;
    }

private:
    static constexpr std::size_t capacity = 256, buckets = 512;
    static constexpr std::uint16_t tombstone = 0xffff;
    struct Segment {
        std::int64_t key = 0;
        std::uint64_t stamp = 0;
        std::size_t bucket = 0;
        Instant<Scale> start{};
        double width = 0;
        std::array<std::array<double, 3>, 4> c{};
    };
    struct Endpoint {
        Instant<Scale> time{};
        State state{};
        bool valid = false;
    };
    Model source_;
    Instant<Scale> first_, last_;
    std::array<Segment, capacity> cache_{};
    // Separate compact hash metadata keeps hits cheap and lets every resident
    // interval use any slot. Set collisions must not turn a 128-date workload
    // into full ELP evaluations when all dates fit the 256-interval capacity.
    std::array<std::uint16_t, buckets> index_{};
    std::size_t tombstones_ = 0;
    std::array<Endpoint, 2> endpoints_{};
    std::uint64_t clock_ = 0;
    std::uint64_t source_evaluations_ = 0;

    static bool less(Instant<Scale> a, Instant<Scale> b) {
        return a.day < b.day || (a.day == b.day && a.frac < b.frac);
    }
    static bool equal(Instant<Scale> a, Instant<Scale> b) {
        return a.day == b.day && a.frac == b.frac;
    }
    static Instant<Scale> normalize(Instant<Scale> t) {
        if (!std::isfinite(t.frac))
            throw std::invalid_argument("Hermite ephemeris: fraction must be finite");
        t = spp::detail::normalized(t);
        if (!(t.frac >= 0 && t.frac < 1))
            throw std::out_of_range("Hermite ephemeris: unnormalizable epoch");
        return t;
    }
    Instant<Scale> checked(Instant<Scale> t) const {
        t = normalize(t);
        if (less(t, first_) || less(last_, t))
            throw std::out_of_range("Hermite ephemeris: epoch outside interpolation domain");
        return t;
    }
    static double difference(Instant<Scale> a, Instant<Scale> b) {
        // This subtraction is only used for endpoints within one day.
        return static_cast<double>(a.day-b.day) + (a.frac-b.frac);
    }
    static double coordinate(const Segment& s, Instant<Scale> t) {
        return difference(t, s.start) / s.width;
    }
    State sample(Instant<Scale> t) {
        for (const auto& n : endpoints_)
            if (n.valid && equal(t,n.time)) return n.state;
        ++source_evaluations_;
        const State s = source_(t);
        for (std::size_t q = 0; q < 3; ++q)
            if (!std::isfinite(s.position_km[q]) || !std::isfinite(s.velocity_km_per_day[q]))
                throw std::invalid_argument("Hermite ephemeris: source returned a non-finite state");
        return s;
    }
    static std::size_t bucket(std::int64_t key) {
        std::uint64_t hash=static_cast<std::uint64_t>(key);
        hash^=hash>>30;hash*=0xbf58476d1ce4e5b9ULL;
        hash^=hash>>27;hash*=0x94d049bb133111ebULL;
        hash^=hash>>31;return static_cast<std::size_t>(hash & (buckets-1));
    }
    void reindex() {
        index_.fill(0);tombstones_=0;
        for(std::size_t i=0;i<cache_.size();++i)if(cache_[i].stamp) {
            auto b=bucket(cache_[i].key);
            while(index_[b])b=(b+1)&(buckets-1);
            cache_[i].bucket=b;index_[b]=static_cast<std::uint16_t>(i+1);
        }
    }
    const Segment& segment(Instant<Scale> t) {
        if (clock_ == std::numeric_limits<std::uint64_t>::max()) clear();
        ++clock_;
        int bin = std::min(static_cast<int>(t.frac*IntervalsPerDay), IntervalsPerDay-1);
        // Multiplication can round across a boundary. Define cells with the
        // same rounded endpoints used for sampling, without extrapolation.
        if (bin > 0 && t.frac < static_cast<double>(bin)/IntervalsPerDay) --bin;
        if (bin+1 < IntervalsPerDay && t.frac >= static_cast<double>(bin+1)/IntervalsPerDay) ++bin;
        const std::int64_t key = t.day*IntervalsPerDay + bin;
        if(tombstones_>buckets/4)reindex();
        auto insertion=buckets;
        auto pos=bucket(key);
        for(std::size_t probe=0;probe<buckets;++probe,pos=(pos+1)&(buckets-1)) {
            const auto slot=index_[pos];
            if(slot==0){if(insertion==buckets)insertion=pos;break;}
            if(slot==tombstone){if(insertion==buckets)insertion=pos;continue;}
            auto& s=cache_[slot-1];
            if(s.key==key){s.stamp=clock_;return s;}
        }
        if(insertion==buckets)throw std::logic_error("Hermite ephemeris: cache index invariant failed");
        Segment* victim=&cache_[0];
        for(auto& s:cache_)if(s.stamp<victim->stamp)victim=&s;

        auto left = Instant<Scale>{t.day, static_cast<double>(bin)/IntervalsPerDay};
        auto right = spp::detail::normalized(Instant<Scale>{t.day, static_cast<double>(bin+1)/IntervalsPerDay});
        if (less(left, first_)) left = first_;
        if (less(last_, right)) right = last_;
        const double width = difference(right, left);
        if (!(width > 0) || !std::isfinite(1.0/width))
            throw std::out_of_range("Hermite ephemeris: interpolation interval is too small");
        const State a = sample(left), b = sample(right);
        Segment built;
        built.key = key; built.stamp = clock_; built.start = left;built.bucket=insertion;
        built.width = width;
        for (std::size_t q = 0; q < 3; ++q) {
            const double delta = b.position_km[q]-a.position_km[q];
            const double va = a.velocity_km_per_day[q]*built.width;
            const double vb = b.velocity_km_per_day[q]*built.width;
            built.c[0][q] = a.position_km[q];
            built.c[1][q] = va;
            built.c[2][q] = 3*delta - 2*va - vb;
            built.c[3][q] = -2*delta + va + vb;
            for (const auto& c : built.c)
                if (!std::isfinite(c[q]))
                    throw std::invalid_argument("Hermite ephemeris: non-finite interpolation coefficient");
        }
        // A throwing provider must not publish a half-filled segment.
        endpoints_ = {Endpoint{left,a,true},Endpoint{right,b,true}};
        if(victim->stamp){index_[victim->bucket]=tombstone;++tombstones_;}
        if(index_[insertion]==tombstone)--tombstones_;
        index_[insertion]=static_cast<std::uint16_t>(victim-cache_.data()+1);
        *victim = built;
        return *victim;
    }
};

}  // namespace spp::model
