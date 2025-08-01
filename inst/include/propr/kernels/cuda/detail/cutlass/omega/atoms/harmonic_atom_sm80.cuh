#pragma once 
#include <cute/tensor.hpp>
#include <propr/utils/constants.h>
#include <propr/utils/preprocessor.cuh>


namespace propr {
    namespace cutlass_atoms { 
        PROPR_FORCE_INLINE 
        PROPR_HOST_DEVICE
        __half2 normalize2(__half2 a, __half2 b) {
            const __half2 two = __float2half2_rn(2.0f);
            const __half2 eps = __float2half2_rn(HLF_EPSILON);
            return __h2div(__hmul2(two, a),__hadd2(__hadd2(a, b), eps));
        }

        struct SM80_16x8x16_F32F16F16F32_TN_HARMONIC_1 {
            using DRegisters = float   [4];
            using ARegisters = uint32_t[4];
            using BRegisters = uint32_t[2];
            using CRegisters = float   [4];

            CUTE_HOST_DEVICE 
            static void
            fma(float          & d0, float          & d1, float & d2, float & d3,
                const uint32_t & a0, const uint32_t & a1, 
                const uint32_t & a2, const uint32_t & a3,
                const uint32_t & b0, const uint32_t & b1,
                float const    & c0, float const    & c1, 
                float const    & c2, float const    & c3) {

                __half2 A0 = *reinterpret_cast<const __half2*>(&a0);
                __half2 A1 = *reinterpret_cast<const __half2*>(&a1);
                __half2 A2 = *reinterpret_cast<const __half2*>(&a2);
                __half2 A3 = *reinterpret_cast<const __half2*>(&a3);
                __half2 B0 = *reinterpret_cast<const __half2*>(&b0);
                __half2 B1 = *reinterpret_cast<const __half2*>(&b1);

                A0 = normalize2(A0, B0);
                A1 = normalize2(A1, B0);
                A2 = normalize2(A2, B1);
                A3 = normalize2(A3, B1);
                
            #if defined(CUTE_ARCH_MMA_SM80_ENABLED)
                asm volatile(
                "mma.sync.aligned.m16n8k16.row.col.f32.f16.f16.f32 "
                "{%0,  %1,  %2,  %3},"
                "{%4,  %5,  %6,  %7},"
                "{%8,  %9},"
                "{%10, %11, %12, %13};\n"
                : "=f"(d0), "=f"(d1), "=f"(d2), "=f"(d3)
                :  "r"(a0),  "r"(a1),  "r"(a2),  "r"(a3),
                    "r"(b0),  "r"(b1),
                    "f"(c0),  "f"(c1),  "f"(c2),  "f"(c3));
            #else
                CUTE_INVALID_CONTROL_PATH("Attempting to use SM80_16x8x16_F32F16F16F32_TN without CUTE_ARCH_MMA_SM80_ENABLED");
            #endif
            }
        };


        struct SM80_16x8x16_F32F16F16F32_TN_HARMONIC_2 {
            using DRegisters = float   [4];
            using ARegisters = uint32_t[4];
            using BRegisters = uint32_t[2];
            using CRegisters = float   [4];

            CUTE_HOST_DEVICE 
            static void
            fma(float          & d0, float          & d1, float         & d2, float         & d3,
                const uint32_t & a0, const uint32_t & a1, 
                const uint32_t & a2, const uint32_t & a3,
                const uint32_t & b0, const uint32_t & b1,
                float const    & c0, float const    & c1, float const   & c2, float const   & c3) {
                __half2 A0 = *reinterpret_cast<const __half2*>(&a0);
                __half2 A1 = *reinterpret_cast<const __half2*>(&a1);
                __half2 A2 = *reinterpret_cast<const __half2*>(&a2);
                __half2 A3 = *reinterpret_cast<const __half2*>(&a3);
                __half2 B0 = *reinterpret_cast<const __half2*>(&b0);
                __half2 B1 = *reinterpret_cast<const __half2*>(&b1);

                A0 = __hmul2(A0, A0);
                A1 = __hmul2(A1, A1);
                A2 = __hmul2(A2, A2);
                A3 = __hmul2(A3, A3);
                B0 = __hmul2(B0, B0);
                B1 = __hmul2(B1, B1);

                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a0)) = A0;
                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a1)) = A1;
                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a2)) = A2;
                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a3)) = A3;
                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&b0)) = B0;
                *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&b1)) = B1;
                #if defined(CUTE_ARCH_MMA_SM80_ENABLED)
                    asm volatile(
                    "mma.sync.aligned.m16n8k16.row.col.f32.f16.f16.f32 "
                    "{%0,  %1,  %2,  %3},"
                    "{%4,  %5,  %6,  %7},"
                    "{%8,  %9},"
                    "{%10, %11, %12, %13};\n"
                    : "=f"(d0), "=f"(d1), "=f"(d2), "=f"(d3)
                    :  "r"(a0),  "r"(a1),  "r"(a2),  "r"(a3),
                        "r"(b0),  "r"(b1),
                        "f"(c0),  "f"(c1),  "f"(c2),  "f"(c3));
                #else
                    CUTE_INVALID_CONTROL_PATH("Attempting to use SM80_16x8x16_F32F16F16F32_TN without CUTE_ARCH_MMA_SM80_ENABLED");
                #endif
            }
        };
    }
}