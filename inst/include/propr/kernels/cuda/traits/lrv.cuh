#pragma once

#include <cub/cub.cuh>

#include <propr/kernels/cuda/traits/common.cuh>

namespace propr {
    namespace cuda {
        namespace traits {
            struct lrv_basic : thread_layout_2d<> {
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };
            
            struct lrv_weighted : thread_layout_2d<> {
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };
            
            struct lrv_alpha : thread_layout_2d<> {
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };

            struct lrv_alpha_weighted : thread_layout_2d<> {
                const static cub::CacheLoadModifier  LoadModifer  = cub::LOAD_CG;
                const static cub::CacheStoreModifier StoreModifer = cub::STORE_CG;
            };
        }
    }
}