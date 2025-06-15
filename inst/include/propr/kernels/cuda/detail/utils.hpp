#pragma once
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

template<typename T>
__device__ __forceinline__
T warp_reduce_sum(T val) {
    val = val +  __shfl_down_sync(0xffffffff, val, 16);
    val = val +  __shfl_down_sync(0xffffffff, val,  8);
    val = val +  __shfl_down_sync(0xffffffff, val,  4);
    val = val +  __shfl_down_sync(0xffffffff, val,  2);
    val = val +  __shfl_down_sync(0xffffffff, val,  1);
    return val;
}


template <typename T, int n>
__device__ __forceinline__
void block_reduce_x(T* sdata) {
    __syncthreads();
    if(threadIdx.x < n/2) sdata[threadIdx.x] += sdata[threadIdx.x + n/2];
    block_reduce_x<n/2>(sdata);
}

template <>
__device__ __forceinline__
void block_reduce_x<0, int>(int* sdata){}
