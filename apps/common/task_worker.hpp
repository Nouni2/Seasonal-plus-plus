#pragma once

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>

namespace spp::app {
struct WorkerCancelled {};
struct WorkerCancellation {
    std::stop_token stop;
    const std::atomic<std::uint64_t>* generation;
    std::uint64_t expected;
    bool operator()()const{return stop.stop_requested() || generation->load()!=expected;}
};

// One coalescing UI worker. Result supplies generation and error members.
// Computation runs without the mutex; only the latest generation is published.
// wake() must remain valid until destruction returns (e.g. before GLFW teardown).
template<class Result>
class TaskWorker {
public:
    using Compute=void(*)(Result&,WorkerCancellation);
private:
    std::mutex mutex_;std::condition_variable cv_;
    std::optional<Result> pending_;std::shared_ptr<const Result> ready_;
    std::atomic<std::uint64_t> generation_{0};
    Compute compute_;void(*wake_)();
    std::jthread thread_;
    static Compute checked(Compute compute) {
        if(!compute)throw std::invalid_argument("UI worker requires a computation");
        return compute;
    }
public:
    explicit TaskWorker(Compute compute,void(*wake)()=nullptr)
        :compute_(checked(compute)),wake_(wake),thread_([this](std::stop_token stop) {
        while(!stop.stop_requested()) {
            std::optional<Result> request;
            {
                std::unique_lock lock(mutex_);
                cv_.wait(lock,[&]{return stop.stop_requested() || pending_.has_value();});
                if(stop.stop_requested())break;
                request=std::move(pending_);pending_.reset();
            }
            auto result=std::make_shared<Result>(std::move(*request));
            try {compute_(*result,{stop,&generation_,result->generation});}
            catch(const WorkerCancelled&) {continue;}
            catch(const std::exception& e) {result->error=e.what();}
            {
                std::lock_guard lock(mutex_);
                if(generation_.load()==result->generation)ready_=std::move(result);
            }
            if(wake_)wake_();
        }
    }) {}
    ~TaskWorker() {
        // The stop predicate must change under the same mutex as wait().
        // Otherwise notification can arrive after a false predicate check
        // but before the waiter sleeps, leaving join() blocked indefinitely.
        {std::lock_guard lock(mutex_);thread_.request_stop();}
        cv_.notify_all();thread_.join();
    }
    TaskWorker(const TaskWorker&)=delete;
    TaskWorker& operator=(const TaskWorker&)=delete;
    std::uint64_t request(Result request) {
        std::lock_guard lock(mutex_);const auto id=++generation_;
        request.generation=id;pending_=std::move(request);cv_.notify_one();return id;
    }
    std::shared_ptr<const Result> ready() {std::lock_guard lock(mutex_);return ready_;}
};
} // namespace spp::app
