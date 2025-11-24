#pragma once

#include <chrono>
#include <string>

#include <propr/utils/host_profiler.hpp>

namespace propr {
    namespace profiler {

    bool host_registration_blocked();

    class host_exclusive_scope_timer {
        public:
            explicit host_exclusive_scope_timer(const char* name);
            ~host_exclusive_scope_timer();

        private:
            bool enabled_;
            std::string name_;
            std::chrono::high_resolution_clock::time_point start_;
    };

        #define PROPR_PROFILE_HOST_EXCLUSIVE_CONCAT(a, b) a##b
        #define PROPR_PROFILE_HOST_EXCLUSIVE_MAKE_NAME(a, b) PROPR_PROFILE_HOST_EXCLUSIVE_CONCAT(a, b)
        #define PROPR_PROFILE_HOST_EXCLUSIVE(NAME) ::propr::profiler::host_exclusive_scope_timer PROPR_PROFILE_HOST_EXCLUSIVE_MAKE_NAME(_propr_profile_scope_, __LINE__)(NAME)
    } 
} 
