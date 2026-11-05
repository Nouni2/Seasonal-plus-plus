#include "spp/model/spk.hpp"
#include "../core/vector.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace spp::model::spk {
namespace {
using namespace spp::detail;
[[noreturn]] void invalid(const char* why){throw std::invalid_argument(std::string("spp SPK: ")+why);}
Instant<TDB> epoch(double seconds) {
    const double d=std::floor(seconds/seconds_per_day);
    if(!std::isfinite(d) || std::abs(d)>1e9)invalid("unusable segment epoch");
    return {jd_j2000_day+static_cast<std::int64_t>(d),(seconds-d*seconds_per_day)/seconds_per_day};
}
double seconds(Instant<TDB> t) {
    // The bounded difference is exact before multiplying; no large Julian-date sum.
    if(!std::isfinite(t.frac))invalid("non-finite TDB fraction");
    if(t.day<jd_j2000_day-1000000000 || t.day>jd_j2000_day+1000000000 || std::abs(t.frac)>1e9)
        throw std::out_of_range("spp SPK: unusable TDB date");
    return static_cast<double>(t.day-jd_j2000_day)*seconds_per_day+t.frac*seconds_per_day;
}
std::size_t integer(double x,std::size_t limit) {
    if(!std::isfinite(x) || x<0 || x>static_cast<double>(limit) || std::floor(x)!=x)
        invalid("invalid directory integer");
    return static_cast<std::size_t>(x);
}
}
struct Kernel::Data {
    struct Directory {std::size_t start,end,size,count;double initial,interval,first,last;};
    std::vector<unsigned char> raw;
    std::vector<Segment> metadata;
    std::vector<Directory> directories;
    bool little=false;
    std::uint64_t bits(std::size_t offset,std::size_t count) const {
        if(offset>raw.size() || count>raw.size()-offset)invalid("truncated kernel");
        // memcpy permits unaligned DAF words and compiles to a single load on
        // native-endian machines. Retain byte assembly for foreign endian data.
        if((little && std::endian::native==std::endian::little) || (!little && std::endian::native==std::endian::big)) {
            if(count==8){std::uint64_t n;std::memcpy(&n,raw.data()+offset,8);return n;}
            if(count==4){std::uint32_t n;std::memcpy(&n,raw.data()+offset,4);return n;}
        }
        std::uint64_t n=0;
        for(std::size_t i=0;i<count;++i)n|=static_cast<std::uint64_t>(raw[offset+i])<<(8*(little?i:count-i-1));
        return n;
    }
    int i32(std::size_t offset)const{return std::bit_cast<std::int32_t>(static_cast<std::uint32_t>(bits(offset,4)));}
    double f64(std::size_t offset)const{return std::bit_cast<double>(bits(offset,8));}
    double word(std::size_t address)const {
        if(!address || address>raw.size()/8)invalid("data address outside kernel");
        return f64(8*(address-1));
    }
    std::string text(std::size_t offset,std::size_t count)const {
        if(offset>raw.size() || count>raw.size()-offset)invalid("truncated kernel string");
        std::string result(reinterpret_cast<const char*>(raw.data()+offset),count);
        while(!result.empty() && (result.back()==' ' || result.back()==0))result.pop_back();
        return result;
    }
    explicit Data(const std::filesystem::path& file) {
        std::ifstream input(file,std::ios::binary|std::ios::ate);
        if(!input)throw std::runtime_error("spp SPK: cannot open "+file.string());
        const auto size=input.tellg();
        if(size<1024 || size%1024!=0 || static_cast<std::uintmax_t>(size)>std::numeric_limits<std::size_t>::max())
            invalid("invalid DAF file size");
        raw.resize(static_cast<std::size_t>(size));input.seekg(0);
        if(!input.read(reinterpret_cast<char*>(raw.data()),size))invalid("cannot read complete kernel");
        if(text(0,8)!="DAF/SPK")invalid("expected DAF/SPK file");
        const auto format=text(88,8);
        if(format!="LTL-IEEE" && format!="BIG-IEEE")invalid("unsupported numeric format");
        little=format=="LTL-IEEE";
        if(i32(8)!=2 || i32(12)!=6)invalid("expected SPK ND=2 NI=6");
        const std::string ftp("FTPSTR:\r:\n:\r\n:\r\0:\x81:\x10\xce:ENDFTP",28);
        if(text(699,28)!=ftp)invalid("damaged FTP integrity string");
        const auto free=i32(84);
        if(free<1 || static_cast<std::size_t>(free)>raw.size()/8+1)invalid("invalid DAF free address");
        std::size_t record=static_cast<std::size_t>(i32(76)),previous=0;
        std::vector<bool> visited(raw.size()/1024+1,false);
        while(record) {
            if(record>=visited.size() || record+1>=visited.size() || visited[record])invalid("invalid/cyclic summary chain");
            visited[record]=true;const auto offset=(record-1)*1024;
            const auto next=integer(f64(offset),visited.size()-1);
            if(integer(f64(offset+8),visited.size()-1)!=previous)invalid("inconsistent summary back link");
            const auto count=integer(f64(offset+16),25);
            for(std::size_t i=0;i<count;++i) {
                const auto at=offset+24+40*i;
                const double first=f64(at),last=f64(at+8);
                if(!std::isfinite(first) || !std::isfinite(last) || first>last)invalid("invalid segment interval");
                Segment info{i32(at+16),i32(at+20),i32(at+24),i32(at+28),epoch(first),epoch(last),text(record*1024+40*i,40)};
                const int a=i32(at+32),b=i32(at+36);
                if(a<1 || b<a || b>=free || static_cast<std::size_t>(b)>raw.size()/8)invalid("invalid segment addresses");
                Directory d{static_cast<std::size_t>(a),static_cast<std::size_t>(b),0,0,0,0,seconds(info.first),seconds(info.last)};
                if(info.type==2 || info.type==3) {
                    if(b-a<4)invalid("short Chebyshev segment");
                    d.initial=word(d.end-3);d.interval=word(d.end-2);
                    d.size=integer(word(d.end-1),d.end-d.start);d.count=integer(word(d.end),d.end-d.start);
                    const std::size_t components=info.type==2?3:6;
                    if(!std::isfinite(d.initial) || !std::isfinite(d.interval) || d.interval<=0 || !d.count
                        || d.size<2+components || (d.size-2)%components || d.size>2+components*1024
                        || d.count>(d.end-d.start-3)/d.size || d.size*d.count+4!=d.end-d.start+1)
                        invalid("invalid Chebyshev directory");
                }
                metadata.push_back(std::move(info));directories.push_back(d);
            }
            previous=record;record=next;
        }
        if(previous!=static_cast<std::size_t>(i32(80)) || metadata.empty())invalid("invalid final summary record");
    }
    std::size_t find(int target,double t) const {
        for(std::size_t i=metadata.size();i--;) {
            const auto& s=metadata[i];
            if(s.target==target && t>=directories[i].first && t<=directories[i].last)return i;
        }
        throw std::out_of_range("spp SPK: no covering segment for body "+std::to_string(target));
    }
    struct EpochCache {
        std::array<std::size_t,64> indices{};
        std::array<State,64> values{};
        std::size_t size=0;
    };
    State evaluate(std::size_t index,double t,EpochCache* cache) const {
        if(!cache)return evaluate(index,t);
        for(std::size_t i=0;i<cache->size;++i)if(cache->indices[i]==index)return cache->values[i];
        const auto state=evaluate(index,t);
        if(cache->size<cache->indices.size()) {
            cache->indices[cache->size]=index;cache->values[cache->size++]=state;
        }
        return state;
    }
    State evaluate(std::size_t index,double t) const {
        const auto& s=metadata[index];const auto& d=directories[index];
        if(s.frame!=1 || (s.type!=2 && s.type!=3))invalid("selected segment needs unsupported frame or SPK type");
        double n=std::floor((t-d.initial)/d.interval);
        if(n==static_cast<double>(d.count))n-=1; // inclusive last endpoint
        if(n<0 || n>=static_cast<double>(d.count))invalid("epoch outside segment records");
        const auto record=d.start+static_cast<std::size_t>(n)*d.size;
        const double mid=word(record),radius=word(record+1);
        if(!std::isfinite(mid) || !std::isfinite(radius) || radius<=0)invalid("invalid Chebyshev record epoch");
        const double x=(t-mid)/radius;
        if(std::abs(x)>1+32*std::numeric_limits<double>::epsilon())invalid("record does not cover requested epoch");
        const std::size_t components=s.type==2?3:6,terms=(d.size-2)/components;
        State result{};
        for(std::size_t c=0;c<components;++c) {
            const auto start=record+2+c*terms;
            double b1=0,b2=0,db1=0,db2=0;
            for(std::size_t k=terms;k-->1;) {
                const double b=word(start+k)+2*x*b1-b2;
                const double db=2*b1+2*x*db1-db2;
                b2=b1;b1=b;db2=db1;db1=db;
            }
            const double value=word(start)+x*b1-b2;
            const double rate=(b1+x*db1-db2)*seconds_per_day/radius;
            if(!std::isfinite(value) || !std::isfinite(rate))invalid("non-finite coefficient/result");
            if(c<3){result.position_km[c]=value;if(s.type==2)result.velocity_km_per_day[c]=rate;}
            else result.velocity_km_per_day[c-3]=value*seconds_per_day;
        }
        return result;
    }
    State relative(int target,int center,double t,EpochCache* cache=nullptr)const {
        if(target==center)return {};
        // Resolve the shared ancestor before summing. Moon/Earth remain near
        // 384000 km instead of subtracting two ~1 AU barycentric states.
        struct Link {int id;std::size_t segment;};
        std::array<Link,64> chain{};std::size_t count=0;int body=target;
        for(;;) {
            if(count==chain.size())invalid("target center chain too deep/cyclic");
            for(std::size_t i=0;i<count;++i)if(chain[i].id==body)invalid("cyclic target centers");
            chain[count++]={body,0};if(body==0)break;
            std::size_t seg;
            try{seg=find(body,t);}catch(const std::out_of_range&){break;}
            chain[count-1].segment=seg;body=metadata[seg].center;
        }
        State result{};body=center;
        std::array<int,64> seen{};std::size_t depth=0;
        for(;;) {
            const auto common=std::find_if(chain.begin(),chain.begin()+count,[&](const Link& v){return v.id==body;});
            if(common!=chain.begin()+count) {
                for(auto i=chain.begin();i!=common;++i) {
                    const auto state=evaluate(i->segment,t,cache);
                    result.position_km=result.position_km+state.position_km;
                    result.velocity_km_per_day=result.velocity_km_per_day+state.velocity_km_per_day;
                }
                return result;
            }
            if(depth==seen.size() || std::find(seen.begin(),seen.begin()+depth,body)!=seen.begin()+depth)
                invalid("cyclic/overlong observer center chain");
            seen[depth++]=body;const auto seg=find(body,t);const auto state=evaluate(seg,t,cache);
            result.position_km=result.position_km-state.position_km;
            result.velocity_km_per_day=result.velocity_km_per_day-state.velocity_km_per_day;
            body=metadata[seg].center;
        }
    }
};
Kernel::Kernel(const std::filesystem::path& file):data_(std::make_shared<Data>(file)){}
State Kernel::state(int target,int center,Instant<TDB> t)const{return data_->relative(target,center,seconds(t));}
void Kernel::states(std::span<const int> targets,int center,Instant<TDB> t,std::span<State> output)const {
    if(targets.size()!=output.size())invalid("batch output size differs from target count");
    const auto et=seconds(t);Data::EpochCache cache;
    for(std::size_t i=0;i<targets.size();++i)output[i]=data_->relative(targets[i],center,et,&cache);
}
std::span<const Segment> Kernel::segments()const{return data_->metadata;}
std::size_t Kernel::bytes()const{return data_->raw.size();}
} // namespace spp::model::spk
