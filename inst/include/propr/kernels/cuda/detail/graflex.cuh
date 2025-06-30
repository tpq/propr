#pragma once

#include <cuda_runtime.h>
#include <cub/device/device_scan.cuh>

#include <propr/utils.hpp>
#include <propr/kernels/cuda/detail/utils.hpp>


namespace propr {
    namespace detail {
        namespace cuda {
            __global__ 
            void compute_odd_ratio_init(cub::ScanTileState<int> a_tile_state,
                                        cub::ScanTileState<int> b_tile_state,
                                        cub::ScanTileState<int> c_tile_state,
                                        cub::ScanTileState<int> d_tile_state, 
                                        int blocks_in_grid){
                a_tile_state.InitializeStatus(blocks_in_grid);
                b_tile_state.InitializeStatus(blocks_in_grid);
                c_tile_state.InitializeStatus(blocks_in_grid);
                d_tile_state.InitializeStatus(blocks_in_grid);
            }

            __global__ void
            compute_odd_ratio(cub::ScanTileState<int> a_tile_state,
                              cub::ScanTileState<int> b_tile_state,
                              cub::ScanTileState<int> c_tile_state,
                              cub::ScanTileState<int> d_tile_state,
                              int* __restrict__ A, int a_stride,
                              int* __restrict__ G, int g_stride,
                              int n, 
                              int *a_acc, int *b_acc, int *c_acc, int *d_acc) {
                using scan_op_t             = cub::Sum;
                using tile_prefix_op        = cub::TilePrefixCallbackOp<int, scan_op_t, cub::ScanTileState<int>>;
                using device_scan_storage_t = typename tile_prefix_op::TempStorage;

                using block_reduce_t        = cub::BlockReduce<int, 128>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;


                __shared__ block_scan_storage_t  block_partials;
                
                __shared__ device_scan_storage_t a_device_partials;
                __shared__ device_scan_storage_t b_device_partials;
                __shared__ device_scan_storage_t c_device_partials;
                __shared__ device_scan_storage_t d_device_partials;

                const int tid = threadIdx.x;
                const int bid = blockIdx.x;
                int stride = gridDim.x * blockDim.x;
                int total_pairs = (n * (n - 1)) / 2;

                for (int k = blockDim.x * bid + tid; k < total_pairs; k += stride) {
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

                    int a_blk_acc = block_reduce_t(block_partials).Sum((1 - a_val) * (1 - g_val)); __syncthreads();
                    int b_blk_acc = block_reduce_t(block_partials).Sum((1 - a_val) * g_val);       __syncthreads();
                    int c_blk_acc = block_reduce_t(block_partials).Sum(a_val * (1 - g_val));       __syncthreads();
                    int d_blk_acc = block_reduce_t(block_partials).Sum(a_val * g_val);             __syncthreads();

                    scan_op_t scan_op{};
                    tile_prefix_op a_prefix(a_tile_state, a_device_partials, scan_op);
                    tile_prefix_op b_prefix(b_tile_state, b_device_partials, scan_op);
                    tile_prefix_op c_prefix(c_tile_state, c_device_partials, scan_op);
                    tile_prefix_op d_prefix(d_tile_state, d_device_partials, scan_op);

                    const int tile_idx = a_prefix.GetTileIdx();
                    if(gridDim.x > 1){
                        if (tile_idx == 0) {
                            if (tid == 0                ) a_tile_state.SetInclusive(tile_idx, a_blk_acc);
                            if (tid == 1*PROPR_WARP_SIZE) b_tile_state.SetInclusive(tile_idx, b_blk_acc);
                            if (tid == 2*PROPR_WARP_SIZE) c_tile_state.SetInclusive(tile_idx, c_blk_acc);
                            if (tid == 3*PROPR_WARP_SIZE) d_tile_state.SetInclusive(tile_idx, d_blk_acc);
                        } else {
                            const int warp_id = tid / PROPR_WARP_SIZE;
                            if (warp_id == 0) {
                                int a_exclusive_prefix = a_prefix(a_blk_acc);
                                if (tid == 0) *a_acc = scan_op(a_exclusive_prefix, a_blk_acc);
                            } else if (warp_id == 1){
                                int b_exclusive_prefix = b_prefix(b_blk_acc);
                                if (tid == PROPR_WARP_SIZE) *b_acc = scan_op(b_exclusive_prefix, b_blk_acc);
                            } else if (warp_id == 2){
                                int c_exclusive_prefix = c_prefix(c_blk_acc);
                                if (tid == 2*PROPR_WARP_SIZE) *c_acc = scan_op(c_exclusive_prefix, c_blk_acc);
                            } else if (warp_id == 3){
                                int d_exclusive_prefix = d_prefix(d_blk_acc);
                                if (tid == 3*PROPR_WARP_SIZE) *d_acc = scan_op(d_exclusive_prefix, d_blk_acc);
                            }
                        }
                    } else {
                        if (tid == 0) {
                            *a_acc = a_blk_acc; *b_acc = b_blk_acc;
                            *c_acc = c_blk_acc; *d_acc = d_blk_acc;
                        }
                    }
                }
            }
        }
    }
}