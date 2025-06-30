#pragma once
#include <cuda_runtime.h>
#define PROPR_WARP_SIZE 32U

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


template<typename T>
__device__ inline T block_reduce_sum(T val, bool final_sync=false) {
    __shared__ T shared_val[PROPR_WARP_SIZE];
    const int lane_id   = threadIdx.x % PROPR_WARP_SIZE;
    const int warp_id   = threadIdx.x / PROPR_WARP_SIZE;
    const int num_warps = blockDim.x  / PROPR_WARP_SIZE;

    T warp_val = warp_reduce_sum(val);
    if (lane_id == 0) { shared_val[warp_id] = warp_val; }
    __syncthreads();
    warp_val = (lane_id < num_warps) ? shared_val[lane_id] : T(0);
    T block_val = warp_reduce_sum(warp_val);
    if (final_sync)  __syncthreads();
    return block_val;
}