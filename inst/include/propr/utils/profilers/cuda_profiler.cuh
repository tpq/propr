#pragma once

#include <cuda_runtime.h>

#include <atomic>
#include <string>
#include <utility>
#include <vector>

namespace propr {
    namespace profiler {

    struct record {
        std::string name;
        double milliseconds;
    };

    void set_enabled(bool on);
    bool is_enabled();

    std::vector<record> consume_records();

    class cuda_scope_timer {
    public:
        cuda_scope_timer(const char* name, cudaStream_t stream);
        ~cuda_scope_timer();

    private:
        bool enabled_;
        std::string name_;
        cudaStream_t stream_;
        cudaEvent_t start_{};
        cudaEvent_t stop_{};
    };

        #define PROPR_PROFILE_CUDA_CONCAT(a, b) a##b
        #define PROPR_PROFILE_CUDA_MAKE_NAME(a, b) PROPR_PROFILE_CUDA_CONCAT(a, b)
        #define PROPR_PROFILE_CUDA(NAME, STREAM) ::propr::profiler::cuda_scope_timer PROPR_PROFILE_CUDA_MAKE_NAME(_propr_profile_scope_, __LINE__)(NAME, STREAM)

    }
} 
