#pragma once

#include <type_traits>

#include <cuda_fp16.h>
#include <cuda_bf16.h>

namespace propr {
        template <typename T>
        struct unsigned_type {
            using type = std::conditional_t<sizeof(T)==1, uint8_t,
                         std::conditional_t<sizeof(T)==2, uint16_t,
                         std::conditional_t<sizeof(T)==4, uint32_t, uint64_t>>>;
        };

        template <typename T> 
        using unsigned_type_t = typename unsigned_type<T>::type;


        template <typename T> 
        constexpr bool is_floating_v  = std::is_floating_point_v<T> || 
                                        std::is_same_v<T,__half>    || 
                                        std::is_same_v<T,__nv_bfloat16>;
    
        template <typename T> 
        constexpr bool is_signed_int_v = std::is_integral_v<T> && std::is_signed_v<T>;

        template <typename> struct vec_scalar_type { using type = void; };
        template <> struct vec_scalar_type<float>  { using type = float;};
        template <> struct vec_scalar_type<float2> { using type = float;};
        template <> struct vec_scalar_type<float3> { using type = float;};
        template <> struct vec_scalar_type<float4> { using type = float;};

        template <typename T, typename U>
        static constexpr bool is_vector_of_v = std::is_same_v<typename vec_scalar_type<T>::type, U>;
}