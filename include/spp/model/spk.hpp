#pragma once

#include "spp/core/state.hpp"
#include "spp/core/time.hpp"
#include <filesystem>
#include <memory>
#include <span>
#include <string>

namespace spp::model::spk {

using State=CartesianState<origin::Relative,frame::Icrf>;
struct Segment {
    int target,center,frame,type;
    Instant<TDB> first,last;
    std::string name;
};

// Immutable, shared, in-memory DAF/SPK reader. No downloads or process-global
// kernel pool. Supports IEEE little/big-endian type 2/3 segments in frame 1.
// NAIF target/center IDs define the runtime origin; states are geometric,
// km and km/TDB day. Last applicable segment wins within this file.
class Kernel {
public:
    explicit Kernel(const std::filesystem::path& file);
    State state(int target,int center,Instant<TDB> t) const;
    // Same states as individual calls, sharing evaluated center segments at
    // this epoch. Scratch storage is local; concurrent readers remain safe.
    void states(std::span<const int> targets,int center,Instant<TDB> t,std::span<State> output) const;
    std::span<const Segment> segments() const;
    std::size_t bytes() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};

} // namespace spp::model::spk
