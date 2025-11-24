#include <Rcpp.h>
#include <propr/utils/host_profiler.hpp>
#include <propr/utils/host_exclusive_profiler.hpp>

#include <atomic>
#include <unordered_map>

namespace propr {
    namespace profiler {

        namespace {
            std::atomic<bool> g_host_enabled{false};
            thread_local std::vector<host_record> g_host_records;

            struct host_timer_state {
                bool enabled{false};
                bool running{false};
                double accumulated_ms{0.0};
                std::chrono::high_resolution_clock::time_point start;
            };

            thread_local std::unordered_map<std::string, host_timer_state> g_host_timers;

            inline double elapsed_ms(const std::chrono::high_resolution_clock::time_point& start,
                                    const std::chrono::high_resolution_clock::time_point& end) {
                return std::chrono::duration<double, std::milli>(end - start).count();
            }
        }

        void set_host_enabled(bool on) { g_host_enabled.store(on, std::memory_order_relaxed); }
        bool host_enabled() { return g_host_enabled.load(std::memory_order_relaxed); }

        std::vector<host_record> consume_host_records() {
            std::vector<host_record> out;
            out.swap(g_host_records);
            return out;
        }

        void record_host_profile(std::string name, double milliseconds) {
            g_host_records.push_back(host_record{std::move(name), milliseconds});
        }

        void start_host_timer(const char* name) {
            if (name == nullptr) return;
            if (host_registration_blocked()) return;
            bool enabled = host_enabled();
            if (!enabled) {
                g_host_timers.erase(name);
                return;
            }
            g_host_timers[name] = host_timer_state{
                true,
                true,
                0.0,
                std::chrono::high_resolution_clock::now(),
            };
        }

        void pause_host_timer(const char* name) {
            if (name == nullptr) return;
            auto it = g_host_timers.find(name);
            if (it == g_host_timers.end()) return;
            auto& timer = it->second;
            if (!timer.enabled || !timer.running) return;
            auto now = std::chrono::high_resolution_clock::now();
            timer.accumulated_ms += elapsed_ms(timer.start, now);
            timer.running = false;
        }

        void resume_host_timer(const char* name) {
            if (name == nullptr) return;
            auto it = g_host_timers.find(name);
            if (it == g_host_timers.end()) return;
            auto& timer = it->second;
            if (!timer.enabled || timer.running) return;
            timer.start   = std::chrono::high_resolution_clock::now();
            timer.running = true;
        }

        void stop_host_timer(const char* name) {
            if (name == nullptr) return;
            auto it = g_host_timers.find(name);
            if (it == g_host_timers.end()) return;
            auto event_name = it->first;
            auto timer      = it->second;
            g_host_timers.erase(it);
            if (!timer.enabled) return;
            if (timer.running) {
                auto now = std::chrono::high_resolution_clock::now();
                timer.accumulated_ms += elapsed_ms(timer.start, now);
            }
            record_host_profile(std::move(event_name), timer.accumulated_ms);
        }

        host_scope_timer::host_scope_timer(const char* name)
            : enabled_(host_enabled() && !host_registration_blocked()),
            name_(name),
            start_(std::chrono::high_resolution_clock::now()) {}

        host_scope_timer::~host_scope_timer() {
            if (!enabled_) return;
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start_).count();
            record_host_profile(std::move(name_), ms);
        }
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
