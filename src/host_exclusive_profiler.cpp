#include <propr/utils/host_exclusive_profiler.hpp>

#include <atomic>

namespace propr {
    namespace profiler {

        namespace {
            thread_local int g_host_exclusive_depth{0};
        }

        bool host_registration_blocked() { return g_host_exclusive_depth > 0; }

        host_exclusive_scope_timer::host_exclusive_scope_timer(const char* name)
            : enabled_(host_enabled() && !host_registration_blocked()),
            name_(name),
            start_(std::chrono::high_resolution_clock::now()) {
            if (enabled_) {
                ++g_host_exclusive_depth;
            }
        }

        host_exclusive_scope_timer::~host_exclusive_scope_timer() {
            if (!enabled_) return;
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start_).count();
            record_host_profile(std::move(name_), ms);
            --g_host_exclusive_depth;
        }

    }
}
