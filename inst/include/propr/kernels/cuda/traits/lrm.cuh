#pragma once

#include <cub/cub.cuh>

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {
            struct lrm_basic {
                using P1_Layout = thread_layout_1d<256>;
                using P2_Layout = thread_layout_2d<16,16>;
            
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };
            
            struct lrm_weighted : thread_layout_2d<>{
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };

            struct lrm_alpha {
                using P1_Layout = thread_layout_1d<256>;
                using P2_Layout = thread_layout_2d<16,16>;

                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };

            struct lrm_alpha_weighted : thread_layout_2d<>{
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };
        }
    }
}
