#pragma once

#include <type_traits>

#include <cub/cub.cuh>

#include <propr/utils/common/preprocessor.cuh>
#include <propr/data/traits.cuh>

namespace propr {
    namespace cuda {
        namespace internal {
            namespace thread {
            
                template <
                    cub::CacheLoadModifier MODIFIER,
                    typename T, 
                    typename U
                >
                PROPR_HOST_DEVICE 
                PROPR_FORCE_INLINE 
                T load(U* ptr) {
                    using Scalar = std::remove_cv_t<U>;
                    static_assert(is_vector_of_v<T, Scalar> || std::is_same_v<T, Scalar>, "U must be the scalar element type of T (e.g., T=float4, U=float) ");
                    #ifdef __CUDA_ARCH__
                        return cub::ThreadLoad<MODIFIER>(reinterpret_cast<const T*>(ptr));
                    #else
                        return *ptr;
                    #endif
                };

                template <
                    cub::CacheStoreModifier MODIFIER, 
                   typename T, 
                   typename U
                >
                PROPR_HOST_DEVICE 
                PROPR_FORCE_INLINE 
                void store(U* ptr, const T& val) {
                    using Scalar = std::remove_cv_t<U>;
                    static_assert(is_vector_of_v<T, Scalar> || std::is_same_v<T, Scalar>, "U must be the scalar element type of T (e.g., T=float4, U=float) ");
                    #ifdef __CUDA_ARCH__
                        return cub::ThreadStore<MODIFIER>(reinterpret_cast<T*>(ptr), val);
                    #else
                        return *ptr;
                    #endif
                };
            }
        }
    }
}
