#pragma once

#include <cuda_runtime.h>

#include <cub/device/device_scan.cuh>

#include <propr/data/types.h>
#include <propr/utils/constants.h>

namespace propr {
    namespace detail {
        namespace cuda {
            
            struct UInt4Sum { // should be moved out
                __device__ __forceinline__ uint4 operator()(const uint4& a, const uint4& b) const {
                    return make_uint4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
                }
            };

            __global__ 
            void compute_odd_ratio_init(cub::ScanTileState<uint4> tile_state,
                                        int blocks_in_grid){
                tile_state.InitializeStatus(blocks_in_grid);
            }

            template <int BLOCK_THREADS>
            __global__ void
            compute_odd_ratio(cub::ScanTileState<uint4> tile_state,
                            unsigned char* __restrict__ A, offset_t a_stride,
                            unsigned char* __restrict__ G, offset_t g_stride,
                            offset_t n, 
                            uint4 *acc) {
                using scan_op_t             = UInt4Sum;
                using tile_prefix_op        = cub::TilePrefixCallbackOp<uint4, scan_op_t, cub::ScanTileState<uint4>>;
                using device_scan_storage_t = typename tile_prefix_op::TempStorage;

                using block_reduce_t        = cub::BlockReduce<uint4, BLOCK_THREADS>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;
                
                __shared__ block_scan_storage_t  block_partials;
                __shared__ device_scan_storage_t device_partials;
                
                const int tid = threadIdx.x;
                const int bid = blockIdx.x;
                
                offset_t stride = gridDim.x * blockDim.x;
                offset_t total_pairs = (n * (n - 1)) / 2;
                
                unsigned int a_acc = 0;
                unsigned int b_acc = 0;
                unsigned int c_acc = 0;
                unsigned int d_acc = 0;

                for (offset_t k = blockDim.x * bid + tid; k < total_pairs; k += stride) {
                    const double disc = 4.0*n*n - 4.0*n - 8.0*k + 1.0;
                    const double s = sqrt(disc);
                    const int a = static_cast<int>(floor(((2.0*n + 1.0) - s) / 2.0));
                    const int S = ((a - 1) * (2*n - a)) / 2;
                    const int b = (a + 1) + (k - S);
                    const int i = a - 1;
                    const int j = b - 1;

                    int a_val = A[i*a_stride + j];
                    int g_val = G[i*g_stride + j];

                    a_acc += (1 - a_val) * (1 - g_val);
                    b_acc += (1 - a_val) * g_val;  
                    c_acc += a_val * (1 - g_val);
                    d_acc += a_val * g_val;
                }
                __syncthreads();
                uint4 values = make_uint4(a_acc,b_acc, c_acc, d_acc);

                scan_op_t scan_op{};
                uint4 blk_acc = block_reduce_t(block_partials).Reduce(values, scan_op);

                tile_prefix_op prefix(tile_state, device_partials, scan_op);
                const int tile_idx = prefix.GetTileIdx();
                if(gridDim.x > 1){
                    if (tile_idx == 0) {
                        if (tid == 0) tile_state.SetInclusive(tile_idx, blk_acc);
                    } else {
                        const unsigned int warp_id = tid / PROPR_WARP_SIZE;
                        if (warp_id == 0) {
                            uint4 exclusive_prefix = prefix(blk_acc);
                            if (tid == 0) {
                                uint4 value = scan_op(exclusive_prefix, blk_acc);
                                if (bid == gridDim.x - 1) *acc = value;
                            }
                        }                  
                    }
                } else {
                    if (tid == 0) *acc = blk_acc;
                }
            }
        }
    }
}