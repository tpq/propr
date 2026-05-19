#pragma once

#include <cub/cub.cuh>
#include <cuda_runtime.h>
#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
            namespace traits {


            struct pairwise_config {
                inline static constexpr int BLK_M = 128;
                inline static constexpr int BLK_K = 8;
                inline static constexpr int TH_Y  = 8;
                inline static constexpr int TH_X  = 8;

                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };

            struct centerNumericMatrix_config: thread_layout_1d<> {};
            struct coordToIndex_config       : thread_layout_1d<> {};
            struct alrRcpp_config            : thread_layout_1d<> {};
            struct indexToCoord_config       : thread_layout_1d<> {};
            struct lltRcpp_config            : thread_layout_1d<> {};
            struct urtRcpp_config            : thread_layout_1d<> {};
            struct labRcpp_config            : thread_layout_1d<> {};
            struct half2mat_config           : thread_layout_1d<> {};
            struct vector2mat_config         : thread_layout_1d<> {};
            struct ratiosRcpp_config         : thread_layout_1d<> {};

            struct cov_config : pairwise_config {};
            struct cor_config : pairwise_config {};
            struct lin_config : pairwise_config {};
            struct vlr_config : pairwise_config {};
            struct phi_config : pairwise_config {};
            struct rho_config : pairwise_config {};

            struct sym_config {
                static constexpr int TILE  = 32; // tile has to be integral multiple of BLK_N
                static constexpr int BLK_N = 16;
            };

        }
    }
}