#pragma once

#include <cuda_runtime.h>

#include <cub/cub.cuh>

#include <propr/utils/constants.h>
#include <propr/utils/preprocessor.cuh>


namespace propr {
    namespace detail {
        namespace cuda {

            template<int BLK_X>
            __global__
            void wtm(float * out,
                     float * __restrict__ x, 
                     float * __restrict__ w,
                     int n){
                static_assert(IS_POWER_OF_2(BLK_X), "BLK_X must be a power of 2");
                using block_reduce_t = cub::BlockReduce<float2, BLK_X>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;
                
                __shared__ block_scan_storage_t partials;

                struct Float2Sum {
                    __device__ __forceinline__ 
                    float2 operator()(const float2& a, const float2& b) const {
                        return make_float2(a.x + b.x, a.y + b.y);
                    }
                };


                float sum_xw_local = 0; float sum_w_local = 0; 
                PROPR_UNROLL
                for (int i = threadIdx.x; i < n; i += BLK_X) {
                    sum_xw_local += x[i] * w[i];
                    sum_w_local  += w[i];
                }
                float2 result = block_reduce_t(partials).Reduce(make_float2(sum_xw_local, sum_w_local), Float2Sum{}); __syncthreads();
                if (threadIdx.x == 0 ){
                    *out =  result.x / result.y;
                }
            };

            template<int BLK_X>
            __global__
            void wtv(float * out,
                     float * __restrict__ x, 
                     float * __restrict__ w,
                     int n){
                static_assert(IS_POWER_OF_2(BLK_X), "BLK_X must be a power of 2");

                using block_reduce_t        = cub::BlockReduce<float4, BLK_X>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;
                
                __shared__ block_scan_storage_t partials;

                struct Combiner {
                    __device__ __forceinline__ 
                    float4 operator()(const float4& a, const float4& b) const {
                        float W_a = a.x, mean_a = a.y, S_a = a.z, W2_a = a.w;
                        float W_b = b.x, mean_b = b.y, S_b = b.z, W2_b = b.w;

                        float W_total = W_a + W_b;
                        float W2_total = W2_a + W2_b;

                        float mean_total = 0.0f;
                        if (W_total != 0.0f) {
                            mean_total = (W_a * mean_a + W_b * mean_b) / W_total;
                        }
                        float delta = mean_b - mean_a;
                        float S_total = S_a + S_b;
                        if (W_total != 0.0f) {
                            S_total += (delta * delta) * (W_a * W_b) / W_total;
                        }
                        return make_float4(W_total, mean_total, S_total, W2_total);
                    }
                };

                float s_local      = 0; 
                float sum_w_local  = 0;
                float sum_w2_local = 0; 
                float mean_local   = 0;
                
                PROPR_UNROLL
                for (int i = threadIdx.x; i < n; i += BLK_X) {
                    float wi       = w[i];
                    float xi       = x[i];
                    float mean_old = mean_local;

                    sum_w_local  += wi;
                    sum_w2_local += wi * wi;
                    mean_local    = mean_local + (wi / sum_w_local) * (xi - mean_old);
                    s_local       = s_local + wi * (xi - mean_old) * (xi - mean_local);
                }

                float4 result  = block_reduce_t(partials).Reduce(make_float4(sum_w_local,mean_local,s_local,sum_w2_local), Combiner{}); 
                __syncthreads();

                if (threadIdx.x == 0) {
                    float denom = result.x * result.x - result.w;
                    if (denom > 0.0f) {
                        *out = result.z * (result.x / denom);
                    } else {
                        *out = NAN;
                    }
                }
            };
        }
    }
}