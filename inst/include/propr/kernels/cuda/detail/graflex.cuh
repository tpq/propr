#pragma once

#include <cuda_runtime.h>
#include <cub/device/device_scan.cuh>
#include <cooperative_groups.h>
#include <propr/utils.hpp>
#include <propr/kernels/cuda/detail/utils.hpp>


namespace propr {
    namespace detail {
        namespace cuda {
            
            struct Int4Sum { // should be moved out
                __device__ __forceinline__ int4 operator()(const int4& a, const int4& b) const {
                    return make_int4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
                }
            };

            //__global__ 
            //void compute_odd_ratio_init(cub::ScanTileState<int4> tile_state,
            //                            int blocks_in_grid){
            //    tile_state.InitializeStatus(blocks_in_grid);
            //}

            __global__ void
            compute_odd_ratio(cub::ScanTileState<int4> tile_state,
                            int* __restrict__ A, int a_stride,
                            int* __restrict__ G, int g_stride,
                            int n, 
                            int4 *acc) {
                using scan_op_t             = Int4Sum;
                using tile_prefix_op        = cub::TilePrefixCallbackOp<int4, scan_op_t, cub::ScanTileState<int4>>;
                using device_scan_storage_t = typename tile_prefix_op::TempStorage;

                using block_reduce_t        = cub::BlockReduce<int4, 512>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;

                __shared__ block_scan_storage_t  block_partials;
                __shared__ device_scan_storage_t device_partials;

                const int tid = threadIdx.x;
                const int bid = blockIdx.x;
                auto grid = cooperative_groups::this_grid();

                int stride = gridDim.x * blockDim.x;
                int total_pairs = (n * (n - 1)) / 2;

                int init_blocks = cub::DivideAndRoundUp(gridDim.x, blockDim.x);

                for (int k = blockDim.x * bid + tid; k < total_pairs; k += stride) {
                    if (blockDim.x*blockIdx.x + threadIdx.x < init_blocks*blockDim.x) tile_state.InitializeStatus(blockDim.x);

                    tile_state.InitializeStatus(blockDim.x);

                    double discriminant = 4.0 * n * n - 4.0 * n - 8.0 * k + 1.0;
                    double sqrt_val     = sqrt(discriminant);
                    double a_float      = ((2.0 * n + 1.0) - sqrt_val) / 2.0;

                    int a = static_cast<int>(floor(a_float));
                    int S = ((a - 1) * (2 * n - a)) / 2;
                    int b = (a + 1) + (k - S);

                    int j = (b - 1);
                    int i = (a - 1);

                    int a_val = A[i*a_stride + j];
                    int g_val = G[i*g_stride + j];

                    // Pack the four values into int4 for vectorized reduction
                    int4 values = make_int4(
                        (1 - a_val) * (1 - g_val),  // a component
                        (1 - a_val) * g_val,        // b component  
                        a_val * (1 - g_val),        // c component
                        a_val * g_val               // d component
                    );
                    
                    int4 blk_acc = block_reduce_t(block_partials).Reduce(values,Int4Sum());
                    __syncthreads();

                    scan_op_t scan_op{};
                    cooperative_groups::sync(grid);
                    tile_prefix_op prefix(tile_state, device_partials, scan_op);

                    const int tile_idx = prefix.GetTileIdx();
                    if(gridDim.x > 1){
                        if (tile_idx == 0) {
                            if (tid == 0) tile_state.SetInclusive(tile_idx, blk_acc);
                        } else {
                            if (tid == 0) {
                                int4 value = scan_op(prefix(blk_acc), blk_acc);
                                acc->x = value.x;
                                acc->y = value.y;
                                acc->z = value.z;
                                acc->w = value.w;
                            }
                        }
                    } else {
                        if (tid == 0) *acc = blk_acc;
                    }
                }
            }
        }
    }
}