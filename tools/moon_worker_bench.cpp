#include "task_worker.hpp"

#include <chrono>
#include <cstdio>
#include <latch>
#include <stdexcept>
#include <string>

namespace {
using Clock=std::chrono::steady_clock;
struct Control {
    std::latch started{1};
    std::atomic<bool> release{false},cancelled{false};
};
struct Result {
    std::uint64_t generation=0;
    std::string error;
    int input=0,output=0;
    Control* control=nullptr;
};
using Worker=spp::app::TaskWorker<Result>;
void require(bool condition,const char* reason){if(!condition)throw std::runtime_error(reason);}
void compute(Result& result,spp::app::WorkerCancellation cancel) {
    if(result.input<0)throw std::runtime_error("controlled calculation error");
    if(result.control) {
        result.control->started.count_down();
        while(!result.control->release.load()) {
            if(cancel()) {result.control->cancelled=true;throw spp::app::WorkerCancelled{};}
            std::this_thread::yield();
        }
    }
    result.output=result.input*2;
}
std::shared_ptr<const Result> ready(Worker& worker,std::uint64_t generation) {
    const auto limit=Clock::now()+std::chrono::seconds(5);
    while(Clock::now()<limit) {
        const auto result=worker.ready();
        if(result && result->generation==generation)return result;
        std::this_thread::yield();
    }
    throw std::runtime_error("worker failed to publish latest request within 5 seconds");
}
void stale(Result& result,spp::app::WorkerCancellation) {
    // An uninterruptible source call can finish after a request is superseded.
    if(result.control) {
        result.control->started.count_down();
        while(!result.control->release.load())std::this_thread::yield();
    }
    result.output=result.input*2;
}
}
int main() {
    try {
        const auto began=Clock::now();
        for(int i=0;i<2048;++i) {Worker worker(compute);}
        const double idle_ms=std::chrono::duration<double,std::milli>(Clock::now()-began).count();
        for(int i=0;i<128;++i) {
            Control control;
            {
                Worker worker(compute);worker.request({0,{},1,0,&control});
                control.started.wait();
            }
            require(control.cancelled,"destruction did not cancel running work");
        }
        for(int i=0;i<128;++i) {
            Control control;Worker worker(compute);
            worker.request({0,{},1,0,&control});control.started.wait();
            std::uint64_t generation=0;
            for(int input=2;input<=100;++input)generation=worker.request({0,{},input,0,nullptr});
            const auto result=ready(worker,generation);
            require(result->input==100 && result->output==200 && result->error.empty() && control.cancelled,
                "coalescing/cancellation published the wrong calculation");
            const auto error_id=worker.request({0,{},-1,0,nullptr});
            require(ready(worker,error_id)->error=="controlled calculation error","calculation exception was lost");
        }
        {
            Control control;Worker worker(stale);
            worker.request({0,{},1,0,&control});control.started.wait();
            const auto generation=worker.request({0,{},7,0,nullptr});control.release=true;
            const auto result=ready(worker,generation);
            require(result->input==7 && result->output==14,"superseded uninterruptible result replaced latest request");
        }
        std::printf("GUI worker: 2048 idle shutdowns in %.3f ms; 128 active cancellations, 128 coalescing/error cycles and stale-source completion: PASS\n",idle_ms);
        return 0;
    }catch(const std::exception& error){std::fprintf(stderr,"FAIL GUI worker: %s\n",error.what());return 1;}
}
