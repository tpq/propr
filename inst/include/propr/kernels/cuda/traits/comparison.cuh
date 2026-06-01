#pragma once

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {

            template <typename Real>
            struct count_values_beyond_thresholds_config_for : thread_layout_1d<256> {
                using value_t  = Real;
                using cutoff_t = Real;
                using block_count_t = unsigned int;
                using accumulator_t = unsigned long long;

                static constexpr int BLOCKS_PER_SM = 8;
                static constexpr int SHARED_CUTOFF_PAD_INTERVAL = 8; // this was by trail and error
            };
            // keeping this for now as it is the default
            using count_values_beyond_thresholds_config = count_values_beyond_thresholds_config_for<double>;

        }
    }
}
