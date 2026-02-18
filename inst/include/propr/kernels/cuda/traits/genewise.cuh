#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {

            struct genewise_connectivity_stats_config {
                static constexpr int THREADS_PER_BLOCK = 256;
                static constexpr int PAIRS_PER_THREAD  = 4;
            };

            template <typename T>
            struct genewise_theta_stats_config_for : thread_layout_1d<1024> {
                static constexpr int RADIX_BITS = 4;
                static constexpr int RADIX_SIZE = 1 << RADIX_BITS;
                static constexpr int RADIX_MASK = RADIX_SIZE - 1;
                static constexpr int CACHE_CAP_VALUES = 2048;
            };
        }
    }
}
