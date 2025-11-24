#pragma once

#include <chrono>
#include <string>
#include <vector>

namespace propr {
    namespace profiler {

    struct host_record {
        std::string name;
        double milliseconds;
    };

    void set_host_enabled(bool on);
    bool host_enabled();
    std::vector<host_record> consume_host_records();

    class host_scope_timer {
    public:
        explicit host_scope_timer(const char* name);
        ~host_scope_timer();

    private:
        bool enabled_;
        std::string name_;
        std::chrono::high_resolution_clock::time_point start_;
    };

    #define PROPR_PROFILE_HOST_CONCAT(a, b) a##b
    #define PROPR_PROFILE_HOST_MAKE_NAME(a, b) PROPR_PROFILE_HOST_CONCAT(a, b)
    #define PROPR_PROFILE_HOST(NAME) ::propr::profiler::host_scope_timer PROPR_PROFILE_HOST_MAKE_NAME(_propr_profile_scope_, __LINE__)(NAME)

    } 
}
