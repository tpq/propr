#pragma once 
#include <cute/tensor.hpp>
#include <propr/utils/constants.h>
#include <propr/utils/preprocessor.cuh>


__device__ int Ag[128][201] = {0};
__device__ int Bg[128][201] = {0};


CUTE_DEVICE
void printAgBg(){
    if (threadIdx.x == 0){
        for(int thread = 0; thread < 128; thread++) {
            printf("thread %d:\n", thread);                
            for(int start_col = 1; start_col < 201; start_col += 50) {
                int end_col = (start_col + 49 < 201) ? start_col + 49 : 200;                    
                for(int col = start_col; col <= end_col; col++) {
                    printf("%d", Ag[thread][col]);
                    if(col < end_col) printf(", ");
                }
                printf("\n");
            }
            printf("\n");
        }
        printf("\n");
        printf("\n");
        printf("\n");

        for(int thread = 0; thread < 128; thread++) {
            printf("thread %d:\n", thread);                
            for(int start_col = 1; start_col < 201; start_col += 50) {
                int end_col = (start_col + 49 < 201) ? start_col + 49 : 200;                    
                for(int col = start_col; col <= end_col; col++) {
                    printf("%d", Bg[thread][col]);
                    if(col < end_col) printf(", ");
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}


namespace propr {
    namespace cutlass_atoms { 

        PROPR_FORCE_INLINE 
        PROPR_HOST_DEVICE
        __half2 normalize2(__half2 a, __half2 b) {
            const __half2 two       = __float2half2_rn(2.0f);
            const __half2 min_norm  = __float2half2_rn(HLF_MIN);
            __half2 sum = __hadd2(a, b);
            __half2 den = __hmax2(sum, min_norm);
            return __h2div(__hmul2(two, a), den);
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

                Ag[threadIdx.x][int(__half2float(A0.x))] += 1;
                Ag[threadIdx.x][int(__half2float(A0.y))] += 1;

                Ag[threadIdx.x][int(__half2float(A2.x))] += 1;
                Ag[threadIdx.x][int(__half2float(A2.y))] += 1;

                Ag[threadIdx.x][int(__half2float(A1.x))] += 1;
                Ag[threadIdx.x][int(__half2float(A1.y))] += 1;

                Ag[threadIdx.x][int(__half2float(A3.x))] += 1;
                Ag[threadIdx.x][int(__half2float(A3.y))] += 1;

                Bg[threadIdx.x][int(__half2float(B0.x))] += 1;
                Bg[threadIdx.x][int(__half2float(B0.y))] += 1;

                Bg[threadIdx.x][int(__half2float(B1.x))] += 1;
                Bg[threadIdx.x][int(__half2float(B1.y))] += 1;


                // if (threadIdx.x == 4 && blockIdx.x == 0 && blockIdx.y == 0) {
                //     printf(
                //         "                        %.2f   %.2f | %.2f   %.2f \n"
                //         " %.2f  %.2f|%.2f  %.2f|                           \n"
                //         " %.2f  %.2f|%.2f  %.2f|                           \n",
                //         __half2float(B0.x), __half2float(B0.y), __half2float(B1.x), __half2float(B1.y),
                //         __half2float(A0.x), __half2float(A0.y), __half2float(A2.x), __half2float(A2.y),
                //         __half2float(A1.x), __half2float(A1.y), __half2float(A3.x), __half2float(A3.y)
                //     );
                // }


                // A0 = normalize2(A0, B0);
                // A1 = normalize2(A1, B0);
                // A2 = normalize2(A2, B1);
                // A3 = normalize2(A3, B1);

                // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a0)) = A0;
                // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a1)) = A1;
                // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a2)) = A2;
                // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&a3)) = A3;
                // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&b0)) = B0;
            //     // *reinterpret_cast<__half2*>(const_cast<uint32_t*>(&b1)) = B1;

            //     uint32_t A0_bits = reinterpret_cast<const uint32_t&>(A0);
            //     uint32_t A1_bits = reinterpret_cast<const uint32_t&>(A1);
            //     uint32_t A2_bits = reinterpret_cast<const uint32_t&>(A2);
            //     uint32_t A3_bits = reinterpret_cast<const uint32_t&>(A3);
                
            // #if defined(CUTE_ARCH_MMA_SM80_ENABLED)
            //     asm volatile(
            //     "mma.sync.aligned.m16n8k16.row.col.f32.f16.f16.f32 "
            //     "{%0,  %1,  %2,  %3},"
            //     "{%4,  %5,  %6,  %7},"
            //     "{%8,  %9},"
            //     "{%10, %11, %12, %13};\n"
            //     : "=f"(d0), "=f"(d1), "=f"(d2), "=f"(d3)
            //     :  "r"(A0_bits),  "r"(A1_bits),  "r"(A2_bits),  "r"(A3_bits),
            //         "r"(b0),  "r"(b1),
            //         "f"(c0),  "f"(c1),  "f"(c2),  "f"(c3));
            // #else
            //     CUTE_INVALID_CONTROL_PATH("Attempting to use SM80_16x8x16_F32F16F16F32_TN without CUTE_ARCH_MMA_SM80_ENABLED");
            // #endif
                double a[8] = {__half2float(A0.x), __half2float(A1.x),   __half2float(A2.x),  __half2float(A3.x), 
                               __half2float(A0.y), __half2float(A1.y),   __half2float(A2.y),  __half2float(A3.y)};

                double b[4] = {__half2float(B0.x), __half2float(B1.x),  
                               __half2float(B0.y), __half2float(B1.y)};

                double c[4] = {c0, c1, c2, c3};
                double d[4] = {d0, d1, d2, d3};

                double const *rA = reinterpret_cast<double const *>(&a);
                double const *rB = reinterpret_cast<double const *>(&b);
                double const *rC = reinterpret_cast<double const *>(&c);
                double *rD       = reinterpret_cast<double *>(&d);

                  asm volatile(
                    "mma.sync.aligned.m16n8k16.row.col.f64.f64.f64.f64"
                    "{%0,  %1,  %2,  %3},"
                    "{%4,  %5,  %6,  %7,  %8,  %9,  %10, %11},"
                    "{%12, %13, %14, %15},"
                    "{%16, %17, %18, %19};\n"
                    : "=d"(rD[0]), "=d"(rD[1]), "=d"(rD[2]), "=d"(rD[3])
                    :  "d"(rA[0]),  "d"(rA[1]),  "d"(rA[2]),  "d"(rA[3]),
                        "d"(rA[4]),  "d"(rA[5]),  "d"(rA[6]),  "d"(rA[7]),
                        "d"(rB[0]),  "d"(rB[1]),  "d"(rB[2]),  "d"(rB[3]),
                        "d"(rC[0]),  "d"(rC[1]),  "d"(rC[2]),  "d"(rC[3]));

                d0 = d[0];
                d1 = d[1];
                d2 = d[2];
                d3 = d[3]; 
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
                const uint32_t & a0_in, const uint32_t & a1_in, 
                const uint32_t & a2_in, const uint32_t & a3_in,
                const uint32_t & b0_in, const uint32_t & b1_in,
                const float    & c0,    const float    & c1, 
                const float    & c2,    const float    & c3) 
            {
                __half2 A0 = *reinterpret_cast<const __half2*>(&a0_in);
                __half2 A1 = *reinterpret_cast<const __half2*>(&a1_in);
                __half2 A2 = *reinterpret_cast<const __half2*>(&a2_in);
                __half2 A3 = *reinterpret_cast<const __half2*>(&a3_in);
                __half2 B0 = *reinterpret_cast<const __half2*>(&b0_in);
                __half2 B1 = *reinterpret_cast<const __half2*>(&b1_in);

                __half2 A0_pre = A0;
                __half2 A1_pre = A1;
                __half2 A2_pre = A2;
                __half2 A3_pre = A3;
                __half2 B0_pre = B0;
                __half2 B1_pre = B1;

                A0 = normalize2(A0, B0);
                A1 = normalize2(A1, B0);
                A2 = normalize2(A2, B1);
                A3 = normalize2(A3, B1);

                // A0 = __hmul2(__hmul2(A0, A0), B0);
                // A1 = __hmul2(__hmul2(A1, A1), B0);
                // A2 = __hmul2(__hmul2(A2, A2), B1);
                // A3 = __hmul2(__hmul2(A3, A3), B1);


                // uint32_t A0_bits = reinterpret_cast<const uint32_t&>(A0);
                // uint32_t A1_bits = reinterpret_cast<const uint32_t&>(A1);
                // uint32_t A2_bits = reinterpret_cast<const uint32_t&>(A2);
                // uint32_t A3_bits = reinterpret_cast<const uint32_t&>(A3);
                // uint32_t B0_bits = reinterpret_cast<const uint32_t&>(B0);
                // uint32_t B1_bits = reinterpret_cast<const uint32_t&>(B1);


                double a[8] = {__half2float(A0.x)*__half2float(A0.x)*__half2float(B0.x), __half2float(A1.x)*__half2float(A1.x)*__half2float(B0.x),   __half2float(A2.x)*__half2float(A2.x)*__half2float(B1.x),  __half2float(A3.x)*__half2float(A3.x)*__half2float(B1.x), 
                               __half2float(A0.y)*__half2float(A0.y)*__half2float(B0.y), __half2float(A1.y)*__half2float(A1.y)*__half2float(B0.y),   __half2float(A2.y)*__half2float(A2.y)*__half2float(B1.y),  __half2float(A3.y)*__half2float(A3.y)*__half2float(B1.y)};

                double b[4] = {__half2float(B0.x), __half2float(B1.x),  
                               __half2float(B0.y), __half2float(B1.y)};

                 double c[4] = {c0, c1, c2, c3};
                 double d[4] = {d0, d1, d2, d3};

                double const *rA = reinterpret_cast<double const *>(&a);
                double const *rB = reinterpret_cast<double const *>(&b);
                double const *rC = reinterpret_cast<double const *>(&c);
                double *rD       = reinterpret_cast<double *>(&d);

                  asm volatile(
                    "mma.sync.aligned.m16n8k16.row.col.f64.f64.f64.f64"
                    "{%0,  %1,  %2,  %3},"
                    "{%4,  %5,  %6,  %7,  %8,  %9,  %10, %11},"
                    "{%12, %13, %14, %15},"
                    "{%16, %17, %18, %19};\n"
                    : "=d"(rD[0]), "=d"(rD[1]), "=d"(rD[2]), "=d"(rD[3])
                    :  "d"(rA[0]),  "d"(rA[1]),  "d"(rA[2]),  "d"(rA[3]),
                        "d"(rA[4]),  "d"(rA[5]),  "d"(rA[6]),  "d"(rA[7]),
                        "d"(rB[0]),  "d"(rB[1]),  "d"(rB[2]),  "d"(rB[3]),
                        "d"(rC[0]),  "d"(rC[1]),  "d"(rC[2]),  "d"(rC[3]));

                d0 = d[0];
                d1 = d[1];
                d2 = d[2];
                d3 = d[3]; 
                // #if defined(CUTE_ARCH_MMA_SM80_ENABLED)
                //     asm volatile(
                //     "mma.sync.aligned.m16n8k16.row.col.f32.f16.f16.f32 "
                //     "{%0, %1, %2, %3},"
                //     "{%4, %5, %6, %7},"
                //     "{%8, %9},"
                //     "{%10, %11, %12, %13};\n"
                //     : "=f"(d0), "=f"(d1), "=f"(d2), "=f"(d3)
                //     : "r"(A0_bits), "r"(A1_bits), "r"(A2_bits), "r"(A3_bits),
                //     "r"(B0_bits), "r"(B1_bits),
                //     "f"(c0), "f"(c1), "f"(c2), "f"(c3)
                //     );
                // #else
                //     CUTE_INVALID_CONTROL_PATH("Attempting to use SM80_16x8x16_F32F16F16F32_TN without CUTE_ARCH_MMA_SM80_ENABLED");
                // #endif
            }
        };
    }
}