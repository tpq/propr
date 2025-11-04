#pragma once

#include <cub/cub.cuh>
#include <cuda_runtime.h>

namespace propr {
    namespace cuda {
            namespace traits {
            struct common_config {
                inline static constexpr int BLK_M = 128;
                inline static constexpr int BLK_K = 8;
                inline static constexpr int TH_Y  = 8;
                inline static constexpr int TH_X  = 8;

                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };

            struct cov_config : common_config {};
            struct cor_config : common_config {};
            struct lin_config : common_config {};
            struct vlr_config : common_config {};
            struct phi_config : common_config {};
            struct rho_config : common_config {};

            struct sym_config {
                static constexpr int TILE  = 32; // tile has to be integral multiple of BLK_N
                static constexpr int BLK_N = 16;
            };

        }
    }
}