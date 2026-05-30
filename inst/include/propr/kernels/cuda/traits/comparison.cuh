#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {

            struct count_values_beyond_thresholds_config : thread_layout_1d<256> {
                using value_t  = double;
                using cutoff_t = double;
                using block_count_t = unsigned int;
                using accumulator_t = unsigned long long;

                static constexpr int BLOCKS_PER_SM = 8;
                static constexpr int SHARED_CUTOFF_PAD_INTERVAL = 8; // this was by trail and error
            };

        }
    }
}
