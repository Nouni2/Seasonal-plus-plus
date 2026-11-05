#include "spp/model/spk.hpp"
#include "core/vector.hpp"
#include <chrono>
#include <bit>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
using namespace spp::detail;
struct Case{int target,center;double et;Vector position,velocity;};
void require(bool condition,const char* why){if(!condition)throw std::runtime_error(why);}
spp::Instant<spp::TDB> epoch(double et) {
    const double d=std::floor(et/spp::seconds_per_day);
    return {spp::jd_j2000_day+static_cast<std::int64_t>(d),(et-d*spp::seconds_per_day)/spp::seconds_per_day};
}
void synthetic() {
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    const auto path=std::filesystem::temp_directory_path()/("spp-spk-"+std::to_string(stamp)+".bsp");
    struct Cleanup{std::filesystem::path path;~Cleanup(){std::error_code e;std::filesystem::remove(path,e);}} cleanup{path};
    for(bool little:{true,false}) {
        std::vector<unsigned char> raw(4096);
        auto bits=[&](std::size_t at,std::uint64_t value,std::size_t size) {
            for(std::size_t i=0;i<size;++i)raw[at+i]=static_cast<unsigned char>(value>>(8*(little?i:size-i-1)));
        };
        auto text=[&](std::size_t at,std::string_view s){std::copy(s.begin(),s.end(),raw.begin()+at);};
        auto i32=[&](std::size_t at,std::uint32_t x){bits(at,x,4);};
        auto real=[&](std::size_t at,double x){bits(at,std::bit_cast<std::uint64_t>(x),8);};
        auto word=[&](std::size_t address,double x){real(8*(address-1),x);};
        auto save=[&] {
            std::ofstream out(path,std::ios::binary|std::ios::trunc);
            out.write(reinterpret_cast<const char*>(raw.data()),static_cast<std::streamsize>(raw.size()));
            require(out.good(),"cannot write temporary synthetic kernel");
        };
        text(0,"DAF/SPK ");i32(8,2);i32(12,6);i32(76,2);i32(80,2);i32(84,418);
        text(88,little?"LTL-IEEE":"BIG-IEEE");text(699,std::string_view("FTPSTR:\r:\n:\r\n:\r\0:\x81:\x10\xce:ENDFTP",28));
        real(1024,0);real(1032,0);real(1040,2);
        for(int n=0;n<2;++n) {
            const std::size_t at=1048+40*n;
            real(at,0);real(at+8,2);i32(at+16,100+n);i32(at+20,n?100:0);i32(at+24,1);i32(at+28,n?3:2);
            i32(at+32,n?400:385);i32(at+36,n?417:399);
        }
        const std::array<double,15> type2{1,1,10,2,0.5,-3,1,0.25,7,-0.5,1,0,2,11,1};
        const std::array<double,18> type3{1,1,100,3,200,-2,300,1,3,0,-2,0,1,0,0,2,14,1};
        for(std::size_t i=0;i<type2.size();++i)word(385+i,type2[i]);
        for(std::size_t i=0;i<type3.size();++i)word(400+i,type3[i]);
        save();const spp::model::spk::Kernel kernel(path);
        for(double t:{0.0,0.3,1.0,1.9,2.0}) {
            const double u=t-1;
            const auto a=kernel.state(100,0,epoch(t)),b=kernel.state(101,100,epoch(t));
            const Vector p{10+2*u+0.5*(2*u*u-1),-3+u+0.25*(2*u*u-1),7-0.5*u+(2*u*u-1)};
            const Vector v{2+2*u,1+u,-0.5+4*u};
            require(norm(a.position_km-p)<1e-12 && norm((1.0/spp::seconds_per_day)*a.velocity_km_per_day-v)<1e-12,
                "type 2 polynomial/derivative differs");
            require(norm(b.position_km-Vector{100+3*u,200-2*u,300+u})<1e-12
                && norm((1.0/spp::seconds_per_day)*b.velocity_km_per_day-Vector{3,-2,1})<1e-12,"type 3 velocity differs");
            const auto combined=kernel.state(101,0,epoch(t)),inverse=kernel.state(0,101,epoch(t));
            require(norm(combined.position_km-a.position_km-b.position_km)<1e-12
                && norm(combined.position_km+inverse.position_km)<1e-12,"SPK center composition/inverse differs");
            const std::array<int,4> targets{101,100,0,101};std::array<spp::model::spk::State,4> batch;
            kernel.states(targets,100,epoch(t),batch);
            for(std::size_t i=0;i<targets.size();++i) {
                const auto single=kernel.state(targets[i],100,epoch(t));
                require(batch[i].position_km==single.position_km && batch[i].velocity_km_per_day==single.velocity_km_per_day,
                    "SPK batch changed common-center, inverse or repeated-target arithmetic");
            }
        }
        bool rejected_size=false;
        try{const std::array<int,1> targets{100};kernel.states(targets,0,epoch(1),{});}
        catch(const std::invalid_argument&){rejected_size=true;}
        require(rejected_size,"SPK accepted mismatched batch output");
        for(int bad=0;bad<6;++bad) {
            const auto good=raw;
            if(bad==0)i32(8,4);
            if(bad==1)i32(1048+32,0);
            if(bad==2)real(1024,2);
            if(bad==3)word(399,2);
            if(bad==4)raw[700]='X';
            if(bad==5)raw.resize(2048);
            save();bool rejected=false;
            try{const spp::model::spk::Kernel broken(path);}catch(const std::invalid_argument&){rejected=true;}
            require(rejected,"corrupt kernel accepted");raw=good;
        }
        // A segment need not provide its own center's path to the barycenter.
        real(1040,1);std::copy(raw.begin()+1088,raw.begin()+1128,raw.begin()+1048);
        save();const spp::model::spk::Kernel partial(path);
        require(norm(partial.state(101,100,epoch(1)).position_km-Vector{100,200,300})<1e-12,
            "relative segment unnecessarily required an unrelated center");
        const std::array<int,2> targets{101,100};std::array<spp::model::spk::State,2> batch;
        partial.states(targets,100,epoch(1),batch);
        require(batch[0].position_km==Vector{100,200,300} && batch[1].position_km==Vector{},"SPK batch required missing ancestors");
    }
    std::puts("SPK synthetic: both endians, type 2/3, derivatives, centers, endpoints and corrupt input: PASS");
}
}
int main(int argc,char** argv) {
    try {
        require(argc==3,"usage: spp_spk_bench de440s.bsp spk_de440s.csv");
        synthetic();
        const auto start=std::chrono::steady_clock::now();const spp::model::spk::Kernel kernel(argv[1]);
        const double load=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        std::ifstream file(argv[2]);require(file.good(),"cannot open CSPICE states");
        std::vector<Case> cases;double maxp=0,maxv=0,lunar_p=0;bool position_ok=true;
        for(std::string line;std::getline(file,line);) {
            if(line.empty() || line.front()=='#')continue;
            for(char& c:line)if(c==',')c=' ';
            std::istringstream in(line);Case q{};in>>q.target>>q.center>>q.et;
            for(auto& x:q.position)in>>x;
            for(auto& x:q.velocity)in>>x;
            require(!in.fail(),"invalid CSPICE fixture");cases.push_back(q);
            const auto s=kernel.state(q.target,q.center,epoch(q.et));
            maxp=std::max(maxp,norm(s.position_km-q.position));
            if(q.target==301 && q.center==399)lunar_p=std::max(lunar_p,norm(s.position_km-q.position));
            position_ok=position_ok && norm(s.position_km-q.position)<std::max(2e-8,8*std::numeric_limits<double>::epsilon()*norm(q.position));
            maxv=std::max(maxv,norm((1.0/spp::seconds_per_day)*s.velocity_km_per_day-q.velocity));
        }
        require(cases.size()==180,"CSPICE fixture coverage changed");
        for(const auto& q:cases) {
            const std::array<int,4> targets{q.target,q.center,q.target,q.center};std::array<spp::model::spk::State,4> batch;
            kernel.states(targets,q.center,epoch(q.et),batch);
            const auto single=kernel.state(q.target,q.center,epoch(q.et));
            require(batch[0].position_km==single.position_km && batch[0].velocity_km_per_day==single.velocity_km_per_day
                && batch[2].position_km==single.position_km && batch[2].velocity_km_per_day==single.velocity_km_per_day
                && batch[1].position_km==Vector{} && batch[3].velocity_km_per_day==Vector{},"batch fixture state changed");
        }
        std::printf("SPK vs independent CSPICE: 180 states, max position %.9g m, velocity %.9g km/s\n",maxp*1000,maxv);
        std::printf("  Moon/Earth position %.9g m; position gate uses 20 micrometres or 8 epsilon times vector norm.\n",lunar_p*1000);
        require(position_ok && maxv<2e-11,"SPK/CSPICE numerical gate failed");
        std::printf("Kernel load %.3f ms; %zu bytes, %zu segments\n",load,kernel.bytes(),kernel.segments().size());
        double checksum=0;
        const auto begin=std::chrono::steady_clock::now();
        for(int batch=0;batch<128;++batch)for(const auto& q:cases) {
            const auto s=kernel.state(q.target,q.center,epoch(q.et));
            checksum+=s.position_km[0]+s.velocity_km_per_day[1];
        }
        std::printf("SPK states: %.3f us/query, checksum %.12g\n",std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-begin).count()/(128*cases.size()),checksum);
        require(std::isfinite(checksum),"invalid SPK benchmark checksum");
        std::puts("PASS: SPK accuracy/performance gates");return 0;
    }catch(const std::exception& e){std::fprintf(stderr,"FAIL SPK benchmark: %s\n",e.what());return 1;}
}
