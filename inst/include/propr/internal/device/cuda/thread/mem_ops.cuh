#pragma once

#include <type_traits>

#include <cub/cub.cuh>

#include <propr/utils/preprocessor.cuh>



template <typename> struct vec_scalar_type { using type = void; };
template <> struct vec_scalar_type<float>  { using type = float;};
template <> struct vec_scalar_type<float2> { using type = float;};
template <> struct vec_scalar_type<float3> { using type = float;};
template <> struct vec_scalar_type<float4> { using type = float;};

template <typename T, typename U>
static constexpr bool is_vector_of_v = std::is_same_v<typename vec_scalar_type<T>::type, U>;


namespace propr {
    namespace cuda {
        namespace internal {

            template <typename T, typename U, cub::CacheLoadModifier MODIFIER = cub::CacheLoadModifier::LOAD_DEFAULT>
            PROPR_HOST_DEVICE 
            PROPR_FORCE_INLINE 
            T load(U* ptr) {
                static_assert(is_vector_of_v<T, U> || std::is_same_v<T,U>, "U must be the scalar element type of T (e.g., T=float4, U=float) ");
                #ifdef __CUDA_ARCH__
                    return cub::ThreadLoad<MODIFIER>(reinterpret_cast<T*>(ptr));
                #else
                    return *ptr;
                #endif
            };

            template <typename T, typename U, cub::CacheStoreModifier MODIFIER = cub::CacheStoreModifier::STORE_DEFAULT>
            PROPR_HOST_DEVICE 
            PROPR_FORCE_INLINE 
            void store(U* ptr, const T& val) {
                static_assert(is_vector_of_v<T, U> || std::is_same_v<T,U>, "U must be the scalar element type of T (e.g., T=float4, U=float) ");
                #ifdef __CUDA_ARCH__
                    return cub::ThreadStore<MODIFIER>(reinterpret_cast<T*>(ptr), val);
                #else
                    return *ptr;
                #endif
            };
        }
    }
}

