#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <propr/utils.hpp>
#include <propr/kernels/cuda/detail/utils.hpp>
#include <cub/cub.cuh> 


namespace propr {
    namespace detail {
        namespace cuda {
            __global__
            void count_joint_zeros(const int* __restrict__ d_X, int X_stride, int nfeats, int* result) {
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                int j = blockIdx.y * blockDim.y + threadIdx.y;
                if (i < nfeats && j < i) {
                    int pair_index = (i * (i - 1)) / 2 + j;
                    result[pair_index] = d_X[i] + d_X[j];
                }
            }
            
            template<int BLK_X>
            __global__
            void count_per_feature(const float* __restrict__ X, int X_stride, int nsubjs, int nfeats, int* result) {
                static_assert(IS_POWER_OF_2(BLK_X), "BLK_X must be a power of 2");
                using BlockReduce = cub::BlockReduce<int, BLK_X>;
                const int tidx = threadIdx.x;
                const int col  = blockIdx.x;
                const int nsubj_padded = (nsubjs + BLK_X - 1) / BLK_X;
                int z_count       = 0;
                __shared__ typename BlockReduce::TempStorage partials;
                for(int i = 0; i < nsubj_padded; i += BLK_X) {
                    if(tidx + i < nsubjs) {
                        z_count += static_cast<int>(X[(tidx + i)  + col*X_stride] == 0);
                    }
                }
                const int val = BlockReduce(partials).Sum(z_count);
                __syncthreads();
                if(tidx == 0) result[col] = val;
            }
        }
    }
}