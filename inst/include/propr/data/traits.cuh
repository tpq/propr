#pragma once

#include <type_traits>

#include <cuda_fp16.h>
#include <cuda_bf16.h>

#include <propr/utils/common/preprocessor.cuh>

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
        template <> struct vec_scalar_type<double2> { using type = double;};

        template <typename T, typename U>
        static constexpr bool is_vector_of_v = std::is_same_v<typename vec_scalar_type<T>::type, U>;

        template <typename T>
        struct cuda_wide_vector;

        template <>
        struct cuda_wide_vector<float> {
            using type = float4;
            static constexpr int lanes = 4;
        };

        template <>
        struct cuda_wide_vector<double> {
            using type = double2;
            static constexpr int lanes = 2;
        };

        template <typename T>
        using cuda_wide_vector_t = typename cuda_wide_vector<T>::type;

        template <typename T>
        inline constexpr int cuda_wide_lanes_v = cuda_wide_vector<T>::lanes;

        PROPR_HOST_DEVICE 
        PROPR_FORCE_INLINE
        float lane_at(const float4& v, int lane) {
            return reinterpret_cast<const float*>(&v)[lane];
        }

        PROPR_HOST_DEVICE 
        PROPR_FORCE_INLINE
        double lane_at(const double2& v, int lane) {
            return reinterpret_cast<const double*>(&v)[lane];
        }

        PROPR_HOST_DEVICE 
        PROPR_FORCE_INLINE
        float4 make_wide_vec(const float* x) {
            return make_float4(x[0], x[1], x[2], x[3]);
        }

        PROPR_HOST_DEVICE 
        PROPR_FORCE_INLINE
        double2 
        make_wide_vec(const double* x) {
            return make_double2(x[0], x[1]);
        }
}
