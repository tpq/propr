#pragma once

#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/warp/reduction.cuh>

namespace propr {
    namespace cuda {
        namespace internal {
            namespace block {

                template<typename ReduceOp, typename T>
                PROPR_NO_DISCARD
                __device__ PROPR_FORCE_INLINE T block_reduce_dyn(T v, ReduceOp op = ReduceOp()) {
                    __shared__ T warp_vals[32];
                    v = warp::warp_reduce(v, op);
                    int lane = threadIdx.x % 32;
                    int warp = threadIdx.x / 32;
                    if (lane == 0) warp_vals[warp] = v;
                    __syncthreads();
                    T result = ReduceOp::identity;
                    if (warp == 0) {
                        result = (lane < (blockDim.x + 31) / 32) ? warp_vals[lane] : ReduceOp::identity;
                        result = warp::warp_reduce(result, op);
                    }
                    return result;
                }
            } // namespace block
        } // namespace internal
    } // namespace cuda
} // namespace propr