#include <Rcpp.h>
#include <propr/utils/host_profiler.hpp>

#include <atomic>

namespace propr::profiler {

    namespace {
        std::atomic<bool> g_host_enabled{false};
        thread_local std::vector<host_record> g_host_records;
    }

    void set_host_enabled(bool on) { g_host_enabled.store(on, std::memory_order_relaxed); }
    bool host_enabled() { return g_host_enabled.load(std::memory_order_relaxed); }

    std::vector<host_record> consume_host_records() {
        std::vector<host_record> out;
        out.swap(g_host_records);
        return out;
    }

    host_scope_timer::host_scope_timer(const char* name)
        : enabled_(host_enabled()), name_(name), start_(std::chrono::high_resolution_clock::now()) {}

    host_scope_timer::~host_scope_timer() {
        if (!enabled_) return;
        auto end = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start_).count();
        g_host_records.push_back(host_record{std::move(name_), ms});
    }
}

using namespace Rcpp;

// [[Rcpp::export]]
void setHostProfile(bool enable) { propr::profiler::set_host_enabled(enable); }

// [[Rcpp::export]]
DataFrame consumeHostProfile() {
    auto recs = propr::profiler::consume_host_records();
    CharacterVector name(recs.size());
    NumericVector ms(recs.size());
    for (size_t i = 0; i < recs.size(); ++i) {
        name[i] = recs[i].name;
        ms[i]   = recs[i].milliseconds;
    }
    return DataFrame::create(_["name"] = name, _["ms"] = ms, _["stringsAsFactors"] = false);
}
