#include <Rcpp.h>
#include <cuda_runtime.h>

#include <propr/utils/cuda_profiler.cuh>

#include <atomic>
#include <mutex>

namespace propr {
    namespace profiler {

    namespace {
        std::atomic<bool> g_enabled{false};
        thread_local std::vector<record> g_records;
    }

    void set_enabled(bool on) { g_enabled.store(on, std::memory_order_relaxed); }

    bool is_enabled() { return g_enabled.load(std::memory_order_relaxed); }

    std::vector<record> consume_records() {
        std::vector<record> out;
        out.swap(g_records);
        return out;
    }

    cuda_scope_timer::cuda_scope_timer(const char* name, cudaStream_t stream)
        : enabled_(is_enabled()), name_(name), stream_(stream) {
        if (!enabled_) return;
        cudaEventCreate(&start_);
        cudaEventCreate(&stop_);
        cudaEventRecord(start_, stream_);
    }

    cuda_scope_timer::~cuda_scope_timer() {
        if (!enabled_) return;
        cudaEventRecord(stop_, stream_);
        cudaEventSynchronize(stop_);
        float ms = 0.0f;
        cudaEventElapsedTime(&ms, start_, stop_);
        g_records.push_back(record{std::move(name_), static_cast<double>(ms)});
        cudaEventDestroy(start_);
        cudaEventDestroy(stop_);
    }

    } 
} 

using namespace Rcpp;

// [[Rcpp::export]]
void setCudaProfile(bool enable) {
    propr::profiler::set_enabled(enable);
}

// [[Rcpp::export]]
DataFrame consumeCudaProfile() {
    auto recs = propr::profiler::consume_records();
    CharacterVector name(recs.size());
    NumericVector ms(recs.size());
    for (size_t i = 0; i < recs.size(); ++i) {
        name[i] = recs[i].name;
        ms[i]   = recs[i].milliseconds;
    }
    return DataFrame::create(_["name"] = name, _["ms"] = ms, _["stringsAsFactors"] = false);
}
