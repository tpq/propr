#pragma once

#include <cuda_runtime.h>

#include <cub/cub.cuh>

#include <propr/data/types.h>
#include <propr/utils/constants.h>
#include <propr/utils/preprocessor.cuh>


// we almost certainly will have problem with occupancy here 
// as the number of registers used per thread is rather quite high
//
// corRcpp: 120 floats/thread est: 2 blks/sm
// covRcpp: 136 floats/thread est: 2 blks/sm
// vlrRcpp: 168 floats/thread est: 1 blks/sm
//
// it might be worth investigating a producer-consumer warp style approach
// so we donot kill occupancy

namespace propr {
    namespace detail {
        namespace cuda {

            struct cov_config { // TODO: move into trait like system
                const static int BLK_M = 128;
                const static int BLK_K = 8;
                const static int TH_Y  = 8;
                const static int TH_X  = 8;
            };

            struct cor_config { // TODO: move into trait like system
                const static int BLK_M = 128;
                const static int BLK_K = 8;
                const static int TH_Y  = 8;
                const static int TH_X  = 8;
            };


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


            template<int BLK_X, int BLK_Y=1, bool col_major=true>
            __global__
            void col_means(
                     float * __restrict__ out, offset_t out_stride,
                     float * __restrict__   x, offset_t x_stride,
                     size_t rows, size_t cols) 
            {
                if constexpr (!col_major){
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;
                    float mean = 0.0;
                    PROPR_UNROLL
                    for (int r = 0; r < rows; ++r) {
                        mean += (x[r * x_stride + col] - mean) / (r + 1);
                    }
                    out[col * out_stride] = mean;
                } else {
                    constexpr int cols_per_block = BLK_Y;
                    constexpr int warps_x = BLK_X / PROPR_WARP_SIZE;

                    const int tx = threadIdx.x % BLK_X;
                    const int ty = threadIdx.x / BLK_X;

                    const int lane = tx & (PROPR_WARP_SIZE - 1);
                    const int warp_x = tx / PROPR_WARP_SIZE;

                    // shared memory holds one partial per warp per col
                    __shared__ float s_warp_sums[ warps_x * BLK_Y];

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < cols);

                        // each thread accumulates a strided sum down the rows
                        float local = 0.0;
                        if (active_col) {
                            for (int r = tx; r < rows; r += BLK_X) {
                                local += x[r + col * x_stride];
                            }
                        }

                        unsigned mask = 0xFFFFFFFFu;
                        PROPR_UNROLL
                        for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                            local += __shfl_down_sync(mask, local, offset);
                        }

                        if (lane == 0) {
                            s_warp_sums[ty * warps_x + warp_x] = local;
                        }
                        __syncthreads();

                        // warp 0 reduces the warp partials for this column
                        if (warp_x == 0) {
                            float partial = (lane < warps_x) ? s_warp_sums[ty * warps_x + lane] : 0.0;
                            PROPR_UNROLL
                            for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                                partial += __shfl_down_sync(mask, partial, offset);
                            }
                            if (lane == 0 && active_col) {
                                out[col * out_stride] = (float)(partial / (float)rows);
                            }
                        }
                        __syncthreads();
                    }
                }
            };


            template <typename T>
            __global__
            void log_transform(T*   __restrict__ out,
                            const T* __restrict__ X,
                            size_t N) 
            {
                const int tid     = blockIdx.x * blockDim.x  + threadIdx.x;
                const int stride  = blockDim.x * gridDim.x;
                const int chunk_size = (N + stride - 1) / stride;
                PROPR_UNROLL
                for (int k = 0; k < chunk_size; ++k) {
                    const offset_t i = tid + k * stride;
                    out[i] = static_cast<T>((i < N)) * log(X[i]);
                }
            };

            template <typename T>
            __global__
            void log_transform_inplace(T* __restrict__ inout, size_t N) 
            {
                const int tid     = blockIdx.x * blockDim.x  + threadIdx.x;
                const int stride  = blockDim.x * gridDim.x;
                const int chunk_size = (N + stride - 1) / stride;
                PROPR_UNROLL
                for (int k = 0; k < chunk_size; ++k) {
                    const offset_t i = tid + k * stride;
                    inout[i] = static_cast<T>((i < N)) * log(inout[i]);
                }
            };

            template<int BLK_X, int BLK_Y = 1, bool col_major = false>
            __global__
            void centerNumericMatrix(
                float* __restrict__ out, offset_t out_stride,
                const float* __restrict__ x, offset_t x_stride,
                size_t rows, size_t cols)
            {
                if constexpr (!col_major) {
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;

                    float mean = 0.0f;
                    for (size_t r = 0; r < rows; ++r) {
                        float v = x[r * x_stride + col];
                        mean += (v - mean) / (float)(r + 1);
                    }
                    for (size_t r = 0; r < rows; ++r) {
                        float v = x[r * x_stride + col];
                        out[r * out_stride + col] = (v - mean);
                    }
                } else {
                    constexpr int cols_per_block = BLK_Y;
                    constexpr int warps_x = BLK_X / PROPR_WARP_SIZE;

                    const int tx = threadIdx.x % BLK_X;
                    const int ty = threadIdx.x / BLK_X;

                    const int lane = tx & (PROPR_WARP_SIZE - 1);
                    const int warp_x = tx / PROPR_WARP_SIZE;

                    __shared__ float s_warp_sums[warps_x * BLK_Y];
                    __shared__ float s_means[BLK_Y];

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < (int)cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < (int)cols);

                        // reduce sum down rows
                        float local = 0.0f;
                        if (active_col) {
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                local += x[r + col * x_stride];
                            }
                        }

                        unsigned mask = 0xFFFFFFFFu;
                        PROPR_UNROLL
                        for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                            local += __shfl_down_sync(mask, local, offset);
                        }

                        if (lane == 0) {
                            s_warp_sums[ty * warps_x + warp_x] = local;
                        }
                        __syncthreads();

                        if (warp_x == 0) {
                            float partial = (lane < warps_x) ? s_warp_sums[ty * warps_x + lane] : 0.0f;
                            PROPR_UNROLL
                            for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                                partial += __shfl_down_sync(mask, partial, offset);
                            }
                            if (lane == 0) {
                                s_means[ty] = active_col ? (partial / (float)rows) : 0.0f;
                            }
                        }
                        __syncthreads();

                        if (active_col) {
                            const float mean = s_means[ty];
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                float v = x[r + col * x_stride];
                                out[r + col * out_stride] = (v - mean);
                            }
                        }
                        __syncthreads();
                    }
                }
            };


            template <class Config>
            __global__ void corRcpp(
                float* __restrict__ out,
                offset_t out_stride,
                const float* __restrict__ x,
                offset_t x_stride,
                int rows,
                int cols
            ) {
                const int M = rows; 
                const int K = cols; 

                const float* A = x;
                const float* B = x;
                      float* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;
                // if (bx > by) return; // (optional) upper triangle only

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];

                float Sa[Config::TH_Y] = {0.0f};   // running S for A lanes (sum of squared deviations)
                float Sb[Config::TH_X] = {0.0f};   // running S for B lanes
                float mu_a[Config::TH_Y] = {0.0f}; // running mean for A lanes
                float mu_b[Config::TH_X] = {0.0f}; // running mean for B lanes
                float accum[Config::TH_Y][Config::TH_X] = {0.0f}; // running cross term

                float frag_a[2][Config::TH_Y];
                float frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);
                float ldg_a_reg[4 * ldg_num_a];
                float ldg_b_reg[4 * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * 4;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const float* A_base = &A[(Config::BLK_M * by) * x_stride];
                const float* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                auto ld_or_zero = [](const float* __restrict__ p,
                                    int r, int c, int ld,
                                    int max_r, int max_c) {
                    return (r >= 0 && r < max_r && c >= 0 && c < max_c) ? p[OFFSET(r, c, ld)] : 0.0f;
                };

                auto store_if_in_bounds = [&](int r, int c, float v) {
                    if (r < M && c < M) {
                        C[OFFSET(r, c, out_stride)] = v;
                    }
                };

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    const int l = (i / A_TILE_ROW_STRIDE) * 4;

                    ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                    ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                    ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                    ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);

                    As[0][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                    As[0][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                    As[0][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                    As[0][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    Bs[0][row_k][col_m + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                }
                __syncthreads();

                FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[0][0][a_tile_index]);
                FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[0][0][a_tile_index + 64]);
                FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[0][0][b_tile_index]);
                FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[0][0][b_tile_index + 64]);

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;

                            ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                            ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                            ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                            ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = tile_idx + B_TILE_ROW_START + i;
                            const int col_m = B_TILE_COL;

                            ldg_b_reg[l + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                            ldg_b_reg[l + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                            ldg_b_reg[l + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                            ldg_b_reg[l + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;

                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem = K - tile_base;
                    const int k_tile = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][0]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][4]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index + 64]);

                        FETCH_FLOAT4(frag_b[(j + 1) & 1][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index + 64]);

                        const float n = float(tile_base + (j + 1));

                        // Update mu_b and Sb first (independent of thread_y)
                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[j & 1][thread_x];
                            float db = b - mu_b[thread_x];
                            float mu_b_new = mu_b[thread_x] + db / n;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        // Update mu_a, Sa, and compute cross term
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[j & 1][thread_y];
                            float da = a - mu_a[thread_y];
                            float mu_a_new = mu_a[thread_y] + da / n;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;
                            const int row_m = A_TILE_ROW_START + i;
                            As[write_stage_idx][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                            As[write_stage_idx][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                            As[write_stage_idx][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                            As[write_stage_idx][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = B_TILE_ROW_START + i;
                            FETCH_FLOAT4(Bs[write_stage_idx][row_k][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[l]);
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const float n_tail = float(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[last_buf][thread_x];
                            float db = b - mu_b[thread_x];
                            float mu_b_new = mu_b[thread_x] + db / n_tail;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[last_buf][thread_y];
                            float da = a - mu_a[thread_y];
                            float mu_a_new = mu_a[thread_y] + da / n_tail;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index]);
                        FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index + 64]);
                        FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index]);
                        FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index + 64]);
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const float eps = 1e-20f;
                float invsig_a[Config::TH_Y];
                float invsig_b[Config::TH_X];

                PROPR_UNROLL
                for (int i = 0; i < Config::TH_Y; ++i) {
                    invsig_a[i] = (Sa[i] > eps) ? rsqrtf(Sa[i]) : 0.0f;
                }
                PROPR_UNROLL
                for (int j = 0; j < Config::TH_X; ++j) {
                    invsig_b[j] = (Sb[j] > eps) ? rsqrtf(Sb[j]) : 0.0f;
                }

                // Helper to compute r(i,j) with clamping
                auto corr_val = [&](int i, int j) -> float {
                    float r = accum[i][j] * invsig_a[i] * invsig_b[j];
                    r = fmaxf(-1.0f, fminf(1.0f, r)); // Clamp tiny FP drift outside [-1,1]
                    return r;
                };

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                for (int i = 0; i < 4; ++i) {
                    store_if_in_bounds(r0 + i, c0 + 0, corr_val(i, 0));
                    store_if_in_bounds(r0 + i, c0 + 1, corr_val(i, 1));
                    store_if_in_bounds(r0 + i, c0 + 2, corr_val(i, 2));
                    store_if_in_bounds(r0 + i, c0 + 3, corr_val(i, 3));

                    store_if_in_bounds(r0 + i, c0 + 64 + 0, corr_val(i, 4));
                    store_if_in_bounds(r0 + i, c0 + 64 + 1, corr_val(i, 5));
                    store_if_in_bounds(r0 + i, c0 + 64 + 2, corr_val(i, 6));
                    store_if_in_bounds(r0 + i, c0 + 64 + 3, corr_val(i, 7));

                    store_if_in_bounds(r0 + 64 + i, c0 + 0, corr_val(i + 4, 0));
                    store_if_in_bounds(r0 + 64 + i, c0 + 1, corr_val(i + 4, 1));
                    store_if_in_bounds(r0 + 64 + i, c0 + 2, corr_val(i + 4, 2));
                    store_if_in_bounds(r0 + 64 + i, c0 + 3, corr_val(i + 4, 3));

                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 0, corr_val(i + 4, 4));
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 1, corr_val(i + 4, 5));
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 2, corr_val(i + 4, 6));
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 3, corr_val(i + 4, 7));
                }
            }
            

            template <class Config>
            __global__ void covRcpp(
                const int norm_type,
                float* __restrict__ out,
                offset_t out_stride,
                const float* __restrict__ x,
                offset_t x_stride,
                int rows,
                int cols
            ) {
                const int M = rows;
                const int K = cols;
                const float* A = x;
                const float* B = x;
                float* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;
                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];

                float mu_a[Config::TH_Y] = {0.0f};
                float mu_b[Config::TH_X] = {0.0f};
                float accum[Config::TH_Y][Config::TH_X] = {0.0f};
                float frag_a[2][Config::TH_Y] = {0.0f};
                float frag_b[2][Config::TH_X] = {0.0f};

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);

                float ldg_a_reg[4 * ldg_num_a];
                float ldg_b_reg[4 * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * 4;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const float* A_base = &A[(Config::BLK_M * by) * x_stride];
                const float* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;

                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                auto ld_or_zero = [](const float* __restrict__ p, int r, int c, int ld, int max_r, int max_c) {
                    return (r >= 0 && r < max_r && c >= 0 && c < max_c) ? p[OFFSET(r, c, ld)] : 0.0f;
                };

                auto store_if_in_bounds = [&](int r, int c, float v) {
                    if (r < M && c < M) {
                        C[OFFSET(r, c, out_stride)] = v;
                    }
                };

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    const int l = (i / A_TILE_ROW_STRIDE) * 4;

                    As[0][A_TILE_COL + 0][row_m] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                    As[0][A_TILE_COL + 1][row_m] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                    As[0][A_TILE_COL + 2][row_m] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                    As[0][A_TILE_COL + 3][row_m] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;

                    Bs[0][row_k][col_m + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                }

                __syncthreads();

                FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[0][0][a_tile_index]);
                FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[0][0][a_tile_index + 64]);
                FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[0][0][b_tile_index]);
                FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[0][0][b_tile_index + 64]);

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;

                            ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                            ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                            ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                            ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);
                        }

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = tile_idx + B_TILE_ROW_START + i;
                            const int col_m = B_TILE_COL;

                            ldg_b_reg[l + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                            ldg_b_reg[l + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                            ldg_b_reg[l + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                            ldg_b_reg[l + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem = K - tile_base;
                    const int k_tile = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][0]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][4]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index + 64]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index + 64]);

                        const float n = float(tile_base + (j + 1));

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[j & 1][thread_x];
                            float db = b - mu_b[thread_x];
                            mu_b[thread_x] += db / n;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[j & 1][thread_y];
                            float da = a - mu_a[thread_y];
                            mu_a[thread_y] += da / n;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;
                            const int row_m = A_TILE_ROW_START + i;

                            As[write_stage_idx][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                            As[write_stage_idx][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                            As[write_stage_idx][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                            As[write_stage_idx][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                        }

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = B_TILE_ROW_START + i;

                            FETCH_FLOAT4(Bs[write_stage_idx][row_k][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[l]);
                        }

                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const float n_tail = float(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[last_buf][thread_x];
                            float db = b - mu_b[thread_x];
                            mu_b[thread_x] += db / n_tail;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[last_buf][thread_y];
                            float da = a - mu_a[thread_y];
                            mu_a[thread_y] += da / n_tail;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index]);
                        FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index + 64]);
                        FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index]);
                        FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index + 64]);
                    }

                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const int ddof = (norm_type != 0);
                const float denom = float(max(1, K + ddof - 1));

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                for (int i = 0; i < 4; ++i) {
                    // Top 4 rows
                    store_if_in_bounds(r0 + i, c0 + 0, accum[i][0] / denom);
                    store_if_in_bounds(r0 + i, c0 + 1, accum[i][1] / denom);
                    store_if_in_bounds(r0 + i, c0 + 2, accum[i][2] / denom);
                    store_if_in_bounds(r0 + i, c0 + 3, accum[i][3] / denom);
                    store_if_in_bounds(r0 + i, c0 + 64 + 0, accum[i][4] / denom);
                    store_if_in_bounds(r0 + i, c0 + 64 + 1, accum[i][5] / denom);
                    store_if_in_bounds(r0 + i, c0 + 64 + 2, accum[i][6] / denom);
                    store_if_in_bounds(r0 + i, c0 + 64 + 3, accum[i][7] / denom);

                    // Bottom 4 rows (+64)
                    store_if_in_bounds(r0 + 64 + i, c0 + 0, accum[i + 4][0] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 1, accum[i + 4][1] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 2, accum[i + 4][2] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 3, accum[i + 4][3] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 0, accum[i + 4][4] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 1, accum[i + 4][5] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 2, accum[i + 4][6] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 3, accum[i + 4][7] / denom);
                }
            }

            template<int BLK_X, int BLK_Y = 1, bool col_major = false>
            __global__
            void clrRcpp(
                float* __restrict__ out, offset_t out_stride,
                const float* __restrict__ x, offset_t x_stride,
                size_t rows, size_t cols)
            {
                if constexpr (!col_major) {
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;

                    // mean of log2(x) down rows
                    float mean = 0.0f;
                    for (size_t r = 0; r < rows; ++r) {
                        float v = log2f(x[r * x_stride + col]);
                        mean += (v - mean) / (float)(r + 1);
                    }
                    for (size_t r = 0; r < rows; ++r) {
                        float v = x[r * x_stride + col];
                        out[r * out_stride + col] = (v - mean);
                    }
                } else {
                    constexpr int cols_per_block = BLK_Y;
                    constexpr int warps_x = BLK_X / PROPR_WARP_SIZE;

                    const int tx = threadIdx.x % BLK_X;
                    const int ty = threadIdx.x / BLK_X;
                    const int lane = tx & (PROPR_WARP_SIZE - 1);
                    const int warp_x = tx / PROPR_WARP_SIZE;

                    __shared__ float s_warp_sums[warps_x * BLK_Y];
                    __shared__ float s_means[BLK_Y];

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < (int)cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < (int)cols);

                        // 1) reduction of sum(log2(x)) down rows
                        float local = 0.0f;
                        if (active_col) {
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                local += log2f(x[r + col * x_stride]);
                            }
                        }

                        unsigned mask = 0xFFFFFFFFu;
                        PROPR_UNROLL
                        for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                            local += __shfl_down_sync(mask, local, offset);
                        }

                        if (lane == 0) {
                            s_warp_sums[ty * warps_x + warp_x] = local;
                        }
                        __syncthreads();

                        if (warp_x == 0) {
                            float partial = (lane < warps_x) ? s_warp_sums[ty * warps_x + lane] : 0.0f;
                            PROPR_UNROLL
                            for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                                partial += __shfl_down_sync(mask, partial, offset);
                            }
                            if (lane == 0) {
                                s_means[ty] = active_col ? (partial / (float)rows) : 0.0f;
                            }
                        }
                        __syncthreads();

                        if (active_col) {
                            const float m = s_means[ty];
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                float v = x[r + col * x_stride];
                                out[r + col * out_stride] = (v - m);
                            }
                        }
                        __syncthreads();
                    }
                }
            }

            template<int BLK_X, int BLK_Y = 1, bool col_major = false>
            __global__
            void alrRcpp(
                const int ivar,
                float* __restrict__ out, offset_t out_stride,
                const float* __restrict__ x, offset_t x_stride,
                size_t rows, size_t cols)
            {
                const int ivar0 = ivar - 1;

                if constexpr (!col_major) {
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;

                    for (size_t r = 0; r < rows; ++r) {
                        float num = log2f(x[r * x_stride + col]);
                        float den = log2f(x[r * x_stride + ivar0]);
                        out[r * out_stride + col] = (num - den);
                    }
                } else {
                    // tile cols and stride rows with BLK_X threads
                    constexpr int cols_per_block = BLK_Y;

                    const int tx = threadIdx.x % BLK_X;
                    const int ty = threadIdx.x / BLK_X;

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < (int)cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < (int)cols);

                        if (active_col) {
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                float num = log2f(x[r + col   * x_stride]);
                                float den = log2f(x[r + ivar0 * x_stride]);
                                out[r + col * out_stride] = (num - den);
                            }
                        }
                        __syncthreads();
                    }
                }
            }
            
            __global__
            void symRcpp(){

            };

            template <class Config>
            __global__
            void phiRcpp(const bool sym,
                         float* __restrict__ out, offset_t out_stride,
                         const float* __restrict__   x, offset_t x_stride,
                               float* __restrict__ row_sums,
                               float* __restrict__ mu_sum,
                               int rows, int cols,
                               int *barrier){
            };


            template <class Config>
            __global__ void
            vlrRcpp(float* __restrict__ out, offset_t out_stride,
                    const float* __restrict__ x, offset_t x_stride,
                    int rows, int cols)
            {
                const int M = rows;
                const int K = cols;

                const float* A = x;
                const float* B = x;
                    float* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X; 
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y; 
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];

                float S [Config::TH_Y][Config::TH_X] = {0.0f};
                float mu[Config::TH_Y][Config::TH_X] = {0.0f};
                
                __syncthreads();

                float frag_a[2][Config::TH_Y];
                float frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);
                float ldg_a_reg[4 * ldg_num_a];
                float ldg_b_reg[4 * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * 4;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const float* A_base = &A[(Config::BLK_M * by) * x_stride];
                const float* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;


                auto ld_or_zero = [](const float* __restrict__ p, int r, int c, int ld, int max_r, int max_c) {
                    return (r < max_r && c < max_c) ? __logf(p[OFFSET(r, c, ld)]) : 0.0f;
                };

                auto store_if_in_bounds = [&](int r, int c, float v) {
                    if (r < M && c < M) C[OFFSET(r, c, out_stride)] = v;
                };

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m  = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    const int l      = (i / A_TILE_ROW_STRIDE) * 4;

                    ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                    ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                    ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                    ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);

                    As[0][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                    As[0][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                    As[0][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                    As[0][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    Bs[0][row_k][col_m + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                }
                __syncthreads();

                FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[0][0][a_tile_index]);
                FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[0][0][a_tile_index + 64]);

                FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[0][0][b_tile_index]);
                FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[0][0][b_tile_index + 64]);

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m  = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l      = (i / A_TILE_ROW_STRIDE) * 4;

                            ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                            ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                            ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                            ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l      = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k  = tile_idx + B_TILE_ROW_START + i;
                            const int col_m  = B_TILE_COL;

                            ldg_b_reg[l + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                            ldg_b_reg[l + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                            ldg_b_reg[l + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                            ldg_b_reg[l + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem       = K - tile_base;
                    const int k_tile    = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max     = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][0]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][4]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index + 64]);

                        FETCH_FLOAT4(frag_b[(j + 1) & 1][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index + 64]);

                        const int k_cur = tile_base + (j + 1);

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[j & 1][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[j & 1][thread_x];
                                const float xval = a - b;

                                const float old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / (float)k_cur;
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                    + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l     = (i / A_TILE_ROW_STRIDE) * 4;
                            const int row_m = A_TILE_ROW_START + i;

                            As[write_stage_idx][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                            As[write_stage_idx][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                            As[write_stage_idx][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                            As[write_stage_idx][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l     = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = B_TILE_ROW_START + i;
                            FETCH_FLOAT4(Bs[write_stage_idx][row_k][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[l]);
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const int k_tail  = tile_base + k_tile;
                        const int last_buf = (j_max & 1);
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[last_buf][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[last_buf][thread_x];
                                const float xval = a - b;

                                const float old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / (float)k_tail;
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                    + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index]);
                        FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index + 64]);

                        FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index]);
                        FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index + 64]);
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const float denom = (K > 1 ? (float)(K - 1) : 1.0f);

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                for (int i = 0; i < 4; ++i) {
                    store_if_in_bounds(r0 + i,      c0 + 0,       S[i][0]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 1,       S[i][1]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 2,       S[i][2]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 3,       S[i][3]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 64 + 0,  S[i][4]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 64 + 1,  S[i][5]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 64 + 2,  S[i][6]    / denom);
                    store_if_in_bounds(r0 + i,      c0 + 64 + 3,  S[i][7]    / denom);

                    store_if_in_bounds(r0 + 64 + i, c0 + 0,       S[i + 4][0] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 1,       S[i + 4][1] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 2,       S[i + 4][2] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 3,       S[i + 4][3] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 0,  S[i + 4][4] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 1,  S[i + 4][5] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 2,  S[i + 4][6] / denom);
                    store_if_in_bounds(r0 + 64 + i, c0 + 64 + 3,  S[i + 4][7] / denom);
                }
            }


            __global__
            void rhoRcpp() {
                // rho based on vlr
            };

            __global__
            void indexToCoord(
                const int N,
                int * __restrict__ V,
                int * __restrict__ row,
                int * __restrict__ col,
                size_t len
            ){
                PROPR_UNROLL
                for (offset_t i = blockIdx.x * blockDim.x + threadIdx.x; i < len; i += blockDim.x * gridDim.x) {
                    row[i] = (V[i] - 1) % N + 1;
                    col[i] = (V[i] - 1) / N + 1;
                }
            };

            __global__
            void coordToIndex(
                const int N,
                int * __restrict__ out, 
                int * __restrict__ row,
                int * __restrict__ col, 
                size_t len){
                // TODO: Check f4 perf
                PROPR_UNROLL
                for (offset_t k = blockIdx.x * blockDim.x + threadIdx.x; k < len; k += blockDim.x * gridDim.x) {
                    out[k] = (col[k] - 1) * N + row[k];
                }
            };

            template <class Config>
            __global__
            void linRcpp(float* __restrict__ out, offset_t out_stride,
                         const float* __restrict__ rho, int rho_stride,
                         const float* __restrict__ x, offset_t x_stride,
                         int rows, int cols){
                const int M = rows;      // features (nfeats)
                const int K = cols;      // samples  (N_samples)

                const float* A = x;
                const float* B = x;
                    float* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;
                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];

                float Sa[Config::TH_Y] = {0.0f};
                float Sb[Config::TH_X] = {0.0f};
                float mu_a[Config::TH_Y] = {0.0f};
                float mu_b[Config::TH_X] = {0.0f};
                float accum[Config::TH_Y][Config::TH_X] = {0.0f};

                float frag_a[2][Config::TH_Y];
                float frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);
                float ldg_a_reg[4 * ldg_num_a];
                float ldg_b_reg[4 * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * 4;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const float* A_base = &A[(Config::BLK_M * by) * x_stride];
                const float* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                auto ld_or_zero = [](const float* __restrict__ p,
                                    int r, int c, int ld, int max_r, int max_c) {
                    return (r >= 0 && r < max_r && c >= 0 && c < max_c) ? p[OFFSET(r, c, ld)] : 0.0f;
                };

                auto rho_at = [&](int r, int c) -> float {
                    return ld_or_zero(rho, r, c, rho_stride, M, M);
                };

                auto store_if_in_bounds = [&](int r, int c, float v) {
                    if (r < M && c < M) C[OFFSET(r, c, out_stride)] = v;
                };

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    const int l = (i / A_TILE_ROW_STRIDE) * 4;
                    ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                    ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                    ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                    ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);

                    As[0][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                    As[0][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                    As[0][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                    As[0][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    Bs[0][row_k][col_m + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                    Bs[0][row_k][col_m + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                }
                __syncthreads();

                FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[0][0][a_tile_index]);
                FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[0][0][a_tile_index + 64]);
                FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[0][0][b_tile_index]);
                FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[0][0][b_tile_index + 64]);

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;
                            ldg_a_reg[l + 0] = ld_or_zero(A_base, row_m, base_k + 0, x_stride, M, K);
                            ldg_a_reg[l + 1] = ld_or_zero(A_base, row_m, base_k + 1, x_stride, M, K);
                            ldg_a_reg[l + 2] = ld_or_zero(A_base, row_m, base_k + 2, x_stride, M, K);
                            ldg_a_reg[l + 3] = ld_or_zero(A_base, row_m, base_k + 3, x_stride, M, K);
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = tile_idx + B_TILE_ROW_START + i;
                            const int col_m = B_TILE_COL;
                            ldg_b_reg[l + 0] = ld_or_zero(B_base, col_m + 0, row_k, x_stride, M, K);
                            ldg_b_reg[l + 1] = ld_or_zero(B_base, col_m + 1, row_k, x_stride, M, K);
                            ldg_b_reg[l + 2] = ld_or_zero(B_base, col_m + 2, row_k, x_stride, M, K);
                            ldg_b_reg[l + 3] = ld_or_zero(B_base, col_m + 3, row_k, x_stride, M, K);
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;

                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem = K - tile_base;
                    const int k_tile = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][0]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j + 1) & 1][4]) = FETCH_FLOAT4(As[load_stage_idx][(j + 1)][a_tile_index + 64]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j + 1) & 1][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j + 1)][b_tile_index + 64]);

                        const float n = float(tile_base + (j + 1));

                        // Update mu_b, Sb
                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[j & 1][thread_x];
                            float db = b - mu_b[thread_x];
                            float mu_b_new = mu_b[thread_x] + db / n;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        // Update mu_a, Sa, and cross
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[j & 1][thread_y];
                            float da = a - mu_a[thread_y];
                            float mu_a_new = mu_a[thread_y] + da / n;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l = (i / A_TILE_ROW_STRIDE) * 4;
                            const int row_m = A_TILE_ROW_START + i;
                            As[write_stage_idx][A_TILE_COL + 0][row_m] = ldg_a_reg[l + 0];
                            As[write_stage_idx][A_TILE_COL + 1][row_m] = ldg_a_reg[l + 1];
                            As[write_stage_idx][A_TILE_COL + 2][row_m] = ldg_a_reg[l + 2];
                            As[write_stage_idx][A_TILE_COL + 3][row_m] = ldg_a_reg[l + 3];
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * 4;
                            const int row_k = B_TILE_ROW_START + i;
                            FETCH_FLOAT4(Bs[write_stage_idx][row_k][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[l]);
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const float n_tail = float(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const float b = frag_b[last_buf][thread_x];
                            float db = b - mu_b[thread_x];
                            float mu_b_new = mu_b[thread_x] + db / n_tail;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const float a = frag_a[last_buf][thread_y];
                            float da = a - mu_a[thread_y];
                            float mu_a_new = mu_a[thread_y] + da / n_tail;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const float b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index]);
                        FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[(write_stage_idx ^ 1)][0][a_tile_index + 64]);
                        FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index]);
                        FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[(write_stage_idx ^ 1)][0][b_tile_index + 64]);
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const float eps     = 1e-20f;
                const float eps_rho = 1e-7f; // clamp for atanh

                auto atanh_clamped = [&](float x) -> float {
                    float xc = fmaxf(-1.0f + eps_rho, fminf(1.0f - eps_rho, x));
                    return 0.5f * logf((1.0f + xc) / (1.0f - xc));
                };


                float Creg[Config::TH_Y][Config::TH_X];

                PROPR_UNROLL
                for (int i = 0; i < Config::TH_Y; ++i) {
                    const int rr = (i < 4) ? (Config::BLK_M * by + c_block_row + i) : (Config::BLK_M * by + c_block_row + 64 + (i - 4));
                    const float inva = (Sa[i] > eps) ? rsqrtf(Sa[i]) : 0.0f;

                    PROPR_UNROLL
                    for (int j = 0; j < Config::TH_X; ++j) {
                        const int cc = (j < 4) ? (Config::BLK_M * bx + c_block_col + j) : (Config::BLK_M * bx + c_block_col + 64 + (j - 4));

                        const float invb = (Sb[j] > eps) ? rsqrtf(Sb[j]) : 0.0f;
                        float r = accum[i][j] * inva * invb;
                        r = fmaxf(-1.0f, fminf(1.0f, r));

                        float out_val = 1.0f;
                        const float rho_ij = rho_at(rr, cc);
                        if (rr > cc) {
                            // lower triangle: z = atanh(rho)
                            out_val = atanh_clamped(rho_ij);
                        } else  if (rr < cc) {
                            // upper triangle: variance formula
                            const float r2     = fmaxf(r * r, eps);
                            const float rho2   = rho_ij * rho_ij;
                            const float denom  = fmaxf((1.0f - rho2) * r2 * (K - 2.0f), eps);
                            const float num    = (1.0f - r2) * rho2;
                            out_val = num / denom;
                        }

                        Creg[i][j] = out_val;
                    }
                }

                const int r0 = (Config::BLK_M * by) + c_block_row;
                const int c0 = (Config::BLK_M * bx) + c_block_col;

                PROPR_UNROLL
                for (int i = 0; i < 4; ++i) {
                    store_if_in_bounds(r0 + i,       c0 + 0, Creg[i][0]);
                    store_if_in_bounds(r0 + i,       c0 + 1, Creg[i][1]);
                    store_if_in_bounds(r0 + i,       c0 + 2, Creg[i][2]);
                    store_if_in_bounds(r0 + i,       c0 + 3, Creg[i][3]);

                    store_if_in_bounds(r0 + 64 + i,  c0 + 0, Creg[i + 4][0]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 1, Creg[i + 4][1]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 2, Creg[i + 4][2]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 3, Creg[i + 4][3]);

                    store_if_in_bounds(r0 + i,       c0 + 64 + 0, Creg[i][4]);
                    store_if_in_bounds(r0 + i,       c0 + 64 + 1, Creg[i][5]);
                    store_if_in_bounds(r0 + i,       c0 + 64 + 2, Creg[i][6]);
                    store_if_in_bounds(r0 + i,       c0 + 64 + 3, Creg[i][7]);

                    store_if_in_bounds(r0 + 64 + i,  c0 + 64 + 0, Creg[i + 4][4]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 64 + 1, Creg[i + 4][5]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 64 + 2, Creg[i + 4][6]);
                    store_if_in_bounds(r0 + 64 + i,  c0 + 64 + 3, Creg[i + 4][7]);
                }
            };

            __global__
            void lltRcpp(
                     float * out, size_t n,
                     float * __restrict__ X, 
                     offset_t x_stride){
                // TODO: Check f4 perf
                offset_t total_pairs = (n * (n - 1)) / 2;

                PROPR_UNROLL
                for (offset_t k = blockDim.x * blockIdx.x + threadIdx.x; k < total_pairs; k += gridDim.x * blockDim.x) {
                    const double t = sqrt(1.0 + 8.0 * (double)k);
                    const offset_t i = static_cast<offset_t>(floor((1.0 + t) / 2.0));
                    const offset_t prev = i * (i - 1) / 2; 
                    const offset_t j = static_cast<offset_t>(k - prev);
                    out[k] = X[j * x_stride + i];
                }
            };

            __global__
            void urtRcpp(
                     float * __restrict__ out, size_t n,
                     float * __restrict__   X,  offset_t x_stride){
                // TODO: Check f4 perf
                offset_t total_pairs = (n * (n - 1)) / 2;

                PROPR_UNROLL
                for (offset_t k = blockDim.x * blockIdx.x + threadIdx.x; k < total_pairs; k += gridDim.x * blockDim.x) {
                    const double t = sqrt(1.0 + 8.0 * (double)k);
                    const offset_t i = static_cast<offset_t>(floor((1.0 + t) / 2.0));
                    const offset_t prev = i * (i - 1) / 2; 
                    const offset_t j = static_cast<offset_t>(k - prev);
                    out[k] = X[i * x_stride + j];
                }
            };

            __global__
            void labRcpp(
                     int * __restrict__ partner,
                     int * __restrict__ pair,
                     size_t n){

                offset_t total_pairs = (n * (n - 1)) / 2;
                int gtid  = blockDim.x * blockIdx.x + threadIdx.x;
                PROPR_UNROLL
                for (offset_t k = gtid; k < total_pairs; k += gridDim.x * blockDim.x) {
                    const double t = sqrt(1.0 + 8.0 * (double)k);
                    const offset_t i = static_cast<offset_t>(floor((1.0 + t) / 2.0));
                    const offset_t prev = i * (i - 1) / 2; 
                    const offset_t j = k - prev;

                    partner[k] = i + 1;
                    pair[k]    = j + 1;
                }
            };

            __global__
            void half2mat(float*       __restrict__ out, 
                          offset_t out_stride,
                          const float* __restrict__ X, 
                          size_t n) {
                // TODO: Check f4 perf
                offset_t total_pairs = n * (n - 1) / 2;

                PROPR_UNROLL
                for (offset_t k = blockIdx.x * blockDim.x + threadIdx.x; k < total_pairs; k += blockDim.x * gridDim.x) {
                    double t = sqrt(1.0 + 8.0 * (double)k);
                    offset_t i = static_cast<offset_t>(floor((1.0 + t) / 2.0));
                    offset_t prev = (offset_t)i * (i - 1) / 2;
                    offset_t j = static_cast<offset_t>(k - prev);
                    out[j * out_stride + i] = X[k];
                    out[i * out_stride + j] = X[k];
                }
            };

            __global__
            void vector2mat(
                float*      __restrict__ out, 
                offset_t    out_stride,
                const float* __restrict__ X,
                const int* __restrict__ i_vec,
                const int* __restrict__ j_vec,
                size_t ni
            ){
                // TODO: Check f4 perf
                PROPR_UNROLL
                for (offset_t idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ni; idx += blockDim.x * gridDim.x) {
                    out[i_vec[idx] - 1 + (j_vec[idx]-1) * out_stride ] = X[idx];
                    out[j_vec[idx] - 1 + (i_vec[idx]-1) * out_stride ] = X[idx];
                }
            };

            __global__
            void ratiosRcpp(
                float*      __restrict__ out, 
                offset_t    out_stride,
                const float* __restrict__ X,
                offset_t    X_stride,
                size_t nfeats, size_t nsamps){

                const offset_t total_pairs = nfeats * (nfeats - 1) / 2;
                const offset_t total_elems = total_pairs * nsamps;

                PROPR_UNROLL
                for (offset_t idx = blockIdx.x * blockDim.x + threadIdx.x; idx < total_elems; idx += blockDim.x * gridDim.x) {
                    offset_t k = idx / nsamps;
                    offset_t s = static_cast<offset_t>(idx % nsamps); 
                    double t = sqrt(1.0 + 8.0 * (double)k);
                    offset_t i = static_cast<offset_t>(floor((1.0 + t) / 2.0));
                    offset_t prev = i * (i - 1) / 2;
                    offset_t j = static_cast<offset_t>(k - prev);
                    out[s + k * out_stride] = X[s + i * X_stride] / X[s + j * X_stride];
                }
            };

            __global__
            void results2matRcpp(
                float*        __restrict__ out   , offset_t out_stride,
                const float* __restrict__ results, offset_t results_stride,
                float diagonal, size_t n){
                offset_t total_pairs = n * (n - 1) / 2;
                PROPR_UNROLL
                for (offset_t k = blockIdx.x * blockDim.x + threadIdx.x; k < total_pairs; k += blockDim.x * gridDim.x) {
                    // TODO: investigate presorting the input
                    int row = int(results[k + 0*results_stride]) - 1;
                    int col = int(results[k + 1*results_stride]) - 1;
                    out[row +  col * out_stride] = row == col ? diagonal : results[k +  2*results_stride];
                    out[col +  row * out_stride] = row == col ? diagonal : results[k +  2*results_stride];
                }
            };

        }
    }
}