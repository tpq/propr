#pragma once

#include <cuda_runtime.h>
#include <cooperative_groups.h>

#include <cub/cub.cuh>

#include <propr/data/math.cuh>
#include <propr/data/traits.cuh>
#include <propr/data/types.h>
#include <propr/utils/common/constants.h>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>
#include <propr/internal/device/cuda/thread/indexing.cuh>


using namespace propr::cuda::internal;

// we almost certainly will have problem with occupancy here 
// as the number of registers used per thread is rather quite high
//
// corRcpp: 120 floats/thread est: 2 blks/sm : OUTDATED
// covRcpp: 136 floats/thread est: 2 blks/sm : OUTDATED
// vlrRcpp: 168 floats/thread est: 1 blks/sm : OUTDATED
//
// it might be worth investigating a producer-consumer warp style approach
// so we donot kill occupancy

namespace propr {
    namespace detail {
        namespace cuda {

            // TODO: move to traits object
            template<typename Real, int BLK_X>
            __global__
            //__launch_bounds__(BLK_X, 1, 1)
            void wtm(Real * out,
                     Real * __restrict__ x,
                     Real * __restrict__ w,
                     int n){
                static_assert(IS_POWER_OF_2(BLK_X), "BLK_X must be a power of 2");
                using block_reduce_t       = cub::BlockReduce<Real, BLK_X>;
                using block_scan_storage_t = typename block_reduce_t::TempStorage;

                __shared__ block_scan_storage_t partials_xw;
                __shared__ block_scan_storage_t partials_w;

                Real sum_xw_local = Real(0);
                Real sum_w_local  = Real(0);
                PROPR_UNROLL
                for (int i = threadIdx.x; i < n; i += BLK_X) {
                    sum_xw_local += x[i] * w[i];
                    sum_w_local  += w[i];
                }
                Real sum_xw = block_reduce_t(partials_xw).Sum(sum_xw_local);
                __syncthreads();
                Real sum_w = block_reduce_t(partials_w).Sum(sum_w_local);
                __syncthreads();
                if (threadIdx.x == 0 ){
                    *out =  sum_xw / sum_w;
                }
            };

            // TODO: move to traits object
            template<typename Real, int BLK_X>
            __global__
            void wtv(Real * out,
                     Real * __restrict__ x,
                     Real * __restrict__ w,
                     int n){
                static_assert(IS_POWER_OF_2(BLK_X), "BLK_X must be a power of 2");

                struct State {
                    Real W;
                    Real mean;
                    Real S;
                    Real W2;
                };

                using block_reduce_t        = cub::BlockReduce<State, BLK_X>;
                using block_scan_storage_t  = typename block_reduce_t::TempStorage;
                
                __shared__ block_scan_storage_t partials;

                struct Combiner {
                    __device__ __forceinline__ 
                    State operator()(const State& a, const State& b) const {
                        Real W_a = a.W, mean_a = a.mean, S_a = a.S, W2_a = a.W2;
                        Real W_b = b.W, mean_b = b.mean, S_b = b.S, W2_b = b.W2;

                        Real W_total = W_a + W_b;
                        Real W2_total = W2_a + W2_b;

                        Real mean_total = Real(0);
                        if (W_total != Real(0)) {
                            mean_total = (W_a * mean_a + W_b * mean_b) / W_total;
                        }
                        Real delta = mean_b - mean_a;
                        Real S_total = S_a + S_b;
                        if (W_total != Real(0)) {
                            S_total += (delta * delta) * (W_a * W_b) / W_total;
                        }
                        return State{W_total, mean_total, S_total, W2_total};
                    }
                };

                Real s_local      = Real(0);
                Real sum_w_local  = Real(0);
                Real sum_w2_local = Real(0);
                Real mean_local   = Real(0);
                
                PROPR_UNROLL
                for (int i = threadIdx.x; i < n; i += BLK_X) {
                    Real wi       = w[i];
                    Real xi       = x[i];
                    Real mean_old = mean_local;

                    sum_w_local  += wi;
                    sum_w2_local += wi * wi;
                    mean_local    = mean_local + (wi / sum_w_local) * (xi - mean_old);
                    s_local       = s_local + wi * (xi - mean_old) * (xi - mean_local);
                }

                State result  = block_reduce_t(partials).Reduce(State{sum_w_local, mean_local, s_local, sum_w2_local}, Combiner{});
                __syncthreads();
                if (threadIdx.x == 0) {
                    Real denom = result.W * result.W - result.W2;
                    if (denom > Real(0)) {
                        *out = result.S * (result.W / denom);
                    } else {
                        *out = static_cast<Real>(NAN);
                    }
                }
            };

            // TODO: move to traits object
            template<typename Real, int BLK_X, int BLK_Y=1, bool row_major=false>
            __global__
            //__launch_bounds__(BLK_X * BLK_Y, 1, 1)
            void col_means(
                     Real * __restrict__ out, offset_t out_stride,
                     Real * __restrict__   x, offset_t x_stride,
                     int rows, int cols) 
            {
                if constexpr (row_major){
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;
                    Real mean = Real(0);
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

                    const int lane   = tx % PROPR_WARP_SIZE;
                    const int warp_x = tx / PROPR_WARP_SIZE;

                    // shared memory holds one partial per warp per col
                    __shared__ Real s_warp_sums[ warps_x * BLK_Y];

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < cols);

                        // each thread accumulates a strided sum down the rows
                        Real local = Real(0);
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
                            Real partial = (lane < warps_x) ? s_warp_sums[ty * warps_x + lane] : Real(0);
                            PROPR_UNROLL
                            for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                                partial += __shfl_down_sync(mask, partial, offset);
                            }
                            if (lane == 0 && active_col) {
                                out[col * out_stride] = partial / static_cast<Real>(rows);
                            }
                        }
                        __syncthreads();
                    }
                }
            };


            template <typename T>
            __global__
            void log_transform(      T*   __restrict__ out,
                               const T* __restrict__ X,
                               size_t N) 
            {
                const int tid     = blockIdx.x * blockDim.x  + threadIdx.x;
                const int stride  = blockDim.x * gridDim.x;
                const int chunk_size = (N + stride - 1) / stride;
                PROPR_UNROLL
                for (int k = 0; k < chunk_size; ++k) {
                    const offset_t i = tid + k * stride;
                    out[i] = static_cast<T>((i < N)) * propr::math::log_t(X[i]);
                }
            };

            template <typename T>
            __global__
            void log_transform_inplace(T* __restrict__ inout, size_t N)  {
                const int tid     = blockIdx.x * blockDim.x  + threadIdx.x;
                const int stride  = blockDim.x * gridDim.x;
                const int chunk_size = (N + stride - 1) / stride;
                PROPR_UNROLL
                for (int k = 0; k < chunk_size; ++k) {
                    const offset_t i = tid + k * stride;
                    inout[i] = static_cast<T>((i < N)) * propr::math::log_t(inout[i]);
                }
            };

            template<typename Real, int BLK_X, int BLK_Y = 1, bool row_major = false>
            __global__
            void centerNumericMatrix(
                Real* __restrict__ out, offset_t out_stride,
                const Real* __restrict__ x, offset_t x_stride,
                int rows, int cols)
            {
                if constexpr (row_major) {
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;

                    Real mean = Real(0);
                    for (size_t r = 0; r < rows; ++r) {
                        Real v = x[r * x_stride + col];
                        mean += (v - mean) / static_cast<Real>(r + 1);
                    }
                    for (size_t r = 0; r < rows; ++r) {
                        Real v = x[r * x_stride + col];
                        out[r * out_stride + col] = (v - mean);
                    }
                } else {
                    constexpr int cols_per_block = BLK_Y;
                    constexpr int warps_x = BLK_X / PROPR_WARP_SIZE;

                    const int tx = threadIdx.x % BLK_X;
                    const int ty = threadIdx.x / BLK_X;

                    const int lane = tx & (PROPR_WARP_SIZE - 1);
                    const int warp_x = tx / PROPR_WARP_SIZE;

                    __shared__ Real s_warp_sums[warps_x * BLK_Y];
                    __shared__ Real s_means[BLK_Y];

                    for (int base_col = blockIdx.x * cols_per_block;
                        base_col < (int)cols;
                        base_col += gridDim.x * cols_per_block)
                    {
                        const int col = base_col + ty;
                        const bool active_col = (col < (int)cols);

                        // reduce sum down rows
                        Real local = Real(0);
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
                            Real partial = (lane < warps_x) ? s_warp_sums[ty * warps_x + lane] : Real(0);
                            PROPR_UNROLL
                            for (int offset = PROPR_WARP_SIZE / 2; offset > 0; offset /= 2) {
                                partial += __shfl_down_sync(mask, partial, offset);
                            }
                            if (lane == 0) {
                                s_means[ty] = active_col ? (partial / static_cast<Real>(rows)) : Real(0);
                            }
                        }
                        __syncthreads();

                        if (active_col) {
                            const Real mean = s_means[ty];
                            for (int r = tx; r < (int)rows; r += BLK_X) {
                                Real v = x[r + col * x_stride];
                                out[r + col * out_stride] = (v - mean);
                            }
                        }
                        __syncthreads();
                    }
                }
            };

            template <typename Real, class Config>
            __global__
            void corRcpp(
                Real* __restrict__ out, offset_t out_stride,
                Real* __restrict__ x  ,  offset_t x_stride,
                int rows, int cols
            ) {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "BLK_M % TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "BLK_M % TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;
                const int K = cols;
                #pragma nv_diag_default 177

                Real* A = x;
                Real* B = x;
                Real* C = out;
                
                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;  // 16
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;  // 16
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real Sa[Config::TH_Y] = {Real(0)};
                Real Sb[Config::TH_X] = {Real(0)};
                Real mu_a[Config::TH_Y] = {Real(0)};
                Real mu_b[Config::TH_X] = {Real(0)};
                Real accum[Config::TH_Y][Config::TH_X] = {Real(0)};

                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                const int A_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int B_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int A_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_THREADS_PER_ROW;
                const int B_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_THREADS_PER_ROW;

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                Real* A_base = &A[(Config::BLK_M * by) * x_stride];
                Real* B_base = &B[(Config::BLK_M * bx) * x_stride];

                {
                    const int a_m0 = tid / A_THREADS_PER_ROW;
                    const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                        const int m = a_m0 + i;
                        const int l = (i / A_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[l], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            As[0][a_k + lane][m] = ldg_a_reg[l + lane];
                        }
                    }
                }
                {
                    const int b_n0 = tid / B_THREADS_PER_ROW;
                    const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                        const int n = b_n0 + i;
                        const int l = (i / B_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[l], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            Bs[0][b_k + lane][n] = ldg_b_reg[l + lane];
                        }
                    }
                }
                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }
                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                            const int m = a_m0 + i;
                            const int l = (i / A_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[l], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k + tile_idx, x_stride)]));
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                            const int n = b_n0 + i;
                            const int l = (i / B_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[l], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k + tile_idx, x_stride)]));
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;

                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem = K - tile_base;
                    const int k_tile = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const Real n = static_cast<Real>(tile_base + (j + 1));

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[j & 1][thread_x];
                            Real db = b - mu_b[thread_x];
                            Real mu_b_new = mu_b[thread_x] + db / n;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            Real da = a - mu_a[thread_y];
                            Real mu_a_new = mu_a[thread_y] + da / n;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                            const int m = a_m0 + i;
                            const int l = (i / A_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][a_k + lane][m] = ldg_a_reg[l + lane];
                            }
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                            const int n = b_n0 + i;
                            const int l = (i / B_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                Bs[write_stage_idx][b_k + lane][n] = ldg_b_reg[l + lane];
                            }
                        }

                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const Real n_tail = static_cast<Real>(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[last_buf][thread_x];
                            Real db = b - mu_b[thread_x];
                            Real mu_b_new = mu_b[thread_x] + db / n_tail;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            Real da = a - mu_a[thread_y];
                            Real mu_a_new = mu_a[thread_y] + da / n_tail;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const Real eps = Real(1e-20);
                Real invsig_a[Config::TH_Y];
                Real invsig_b[Config::TH_X];

                PROPR_UNROLL
                for (int i = 0; i < Config::TH_Y; ++i) {
                    invsig_a[i] = (Sa[i] > eps) ? propr::math::rsqrt_t(Sa[i]) : static_cast<Real>(PROPR_R_NA_REAL);
                }
                PROPR_UNROLL
                for (int j = 0; j < Config::TH_X; ++j) {
                    invsig_b[j] = (Sb[j] > eps) ? propr::math::rsqrt_t(Sb[j]) : static_cast<Real>(PROPR_R_NA_REAL);
                }

                auto corr_val = [&](int i, int j) -> Real {
                    return accum[i][j] * invsig_a[i] * invsig_b[j];
                };

                PROPR_UNROLL
                for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                    const int row_offset = (thread_y < Config::TH_Y / 2)
                        ? thread_y
                        : Config::BLK_M / 2 + (thread_y - Config::TH_Y / 2);
                    const int row = Config::BLK_M * by + c_block_row + row_offset;

                    PROPR_UNROLL
                    for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                        const int col_offset = (thread_x < Config::TH_X / 2)
                            ? thread_x
                            : Config::BLK_M / 2 + (thread_x - Config::TH_X / 2);
                        const int col = Config::BLK_M * bx + c_block_col + col_offset;
                        thread::store<Config::StoreModifer, Real>(&C[OFFSET(row, col, out_stride)], corr_val(thread_y, thread_x));
                    }
                }
            }


            template <typename Real, class Config>
            __global__ void covRcpp(
                const int norm_type,
                Real* __restrict__ out, offset_t out_stride,
                Real* __restrict__ x, offset_t x_stride,
                int rows, int cols /* K (true) */
            ) {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "BLK_M % TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "BLK_M % TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;
                const int K = cols;
                #pragma nv_diag_default 177

                Real* A = x;
                Real* B = x;
                Real* C = out;

                const int bx = blockIdx.x, by = blockIdx.y;
                const int tx = threadIdx.x, ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real mu_a[Config::TH_Y] = {Real(0)};
                Real mu_b[Config::TH_X] = {Real(0)};
                Real accum[Config::TH_Y][Config::TH_X] = {Real(0)};
                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int A_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int B_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int A_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_THREADS_PER_ROW;
                const int B_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_THREADS_PER_ROW;

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                Real* A_base = &A[(Config::BLK_M * by) * x_stride];
                Real* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                {
                    const int a_m0 = tid / A_THREADS_PER_ROW;
                    const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                        const int m = a_m0 + i;
                        const int l = (i / A_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[l], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            As[0][a_k + lane][m] = ldg_a_reg[l + lane];
                        }
                    }
                }
                {
                    const int b_n0 = tid / B_THREADS_PER_ROW;
                    const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                        const int n = b_n0 + i;
                        const int l = (i / B_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[l], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            Bs[0][b_k + lane][n] = ldg_b_reg[l + lane];
                        }
                    }
                }

                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }
                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                            const int m = a_m0 + i;
                            const int l = (i / A_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[l], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k + tile_idx, x_stride)]));
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                            const int n = b_n0 + i;
                            const int l = (i / B_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[l], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k + tile_idx, x_stride)]));
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem      = K - tile_base;
                    const int k_tile   = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max    = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const Real n = static_cast<Real>(tile_base + (j + 1));

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[j & 1][thread_x];
                            Real db = b - mu_b[thread_x];
                            mu_b[thread_x] += db / n;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            Real da = a - mu_a[thread_y];
                            mu_a[thread_y] += da / n;
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_ROW_STRIDE) {
                            const int m = a_m0 + i;
                            const int l = (i / A_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][a_k + lane][m] = ldg_a_reg[l + lane];
                            }
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_ROW_STRIDE) {
                            const int n = b_n0 + i;
                            const int l = (i / B_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                Bs[write_stage_idx][b_k + lane][n] = ldg_b_reg[l + lane];
                            }
                        }

                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const Real n_tail = static_cast<Real>(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[last_buf][thread_x];
                            Real db = b - mu_b[thread_x];
                            mu_b[thread_x] += db / n_tail;
                        }
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            Real da = a - mu_a[thread_y];
                            mu_a[thread_y] += da / n_tail;
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }

                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const int ddof = (norm_type != 0);
                const Real denom = static_cast<Real>(max(1, K + ddof - 1));

                PROPR_UNROLL
                for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                    const int row_offset = (thread_y < Config::TH_Y / 2)
                        ? thread_y
                        : Config::BLK_M / 2 + (thread_y - Config::TH_Y / 2);
                    const int row = Config::BLK_M * by + c_block_row + row_offset;

                    PROPR_UNROLL
                    for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                        const int col_offset = (thread_x < Config::TH_X / 2)
                            ? thread_x
                            : Config::BLK_M / 2 + (thread_x - Config::TH_X / 2);
                        const int col = Config::BLK_M * bx + c_block_col + col_offset;
                        thread::store<Config::StoreModifer, Real>(&C[OFFSET(row, col, M)], accum[thread_y][thread_x] / denom);
                    }
                }
            }

            template<typename Real, int BLK_X, int BLK_Y = 1, bool row_major = false>
            __global__
            void clrRcpp(
                      Real* __restrict__ out, offset_t out_stride,
                const Real* __restrict__ x, offset_t x_stride,
                int rows, int cols)
            {
                constexpr int WARPS_X = BLK_X / PROPR_WARP_SIZE;
                const int tx    = threadIdx.x % BLK_X;      // 0 .. BLK_X-1
                const int ty    = threadIdx.x / BLK_X;      // 0 .. BLK_Y-1
                const int lane  = tx & (PROPR_WARP_SIZE - 1);
                const int warpX = tx / PROPR_WARP_SIZE;

                const int row = blockIdx.x * BLK_Y + ty;
                const bool active_row = (row < rows);

                __shared__ Real s_warp_sums[BLK_Y * WARPS_X];
                __shared__ Real s_means[BLK_Y];

                Real local = Real(0);
                if (active_row) {
                    for (int c = tx; c < cols; c += BLK_X) {
                        local += propr::math::log_t(x[row + c * x_stride]);
                    }
                }

                unsigned mask = 0xFFFFFFFFu;
                for (int off = PROPR_WARP_SIZE/2; off > 0; off >>= 1) {
                    local += __shfl_down_sync(mask, local, off);
                }

                if (lane == 0) {
                    s_warp_sums[ty * WARPS_X + warpX] = local;
                }
                __syncthreads();

                if (warpX == 0) {
                    Real partial = (lane < WARPS_X) ? s_warp_sums[ty * WARPS_X + lane] : Real(0);
                    for (int off = PROPR_WARP_SIZE/2; off > 0; off >>= 1) {
                        partial += __shfl_down_sync(mask, partial, off);
                    }
                    if (lane == 0) {
                        s_means[ty] = active_row ? (partial / static_cast<Real>(cols)) : Real(0);
                    }
                }
                __syncthreads();

                if (active_row) {
                    const Real m = s_means[ty];
                    for (int c = tx; c < cols; c += BLK_X) {
                        Real v = propr::math::log_t(x[row + c * x_stride]);
                        out[row + c * out_stride] = v - m;
                    }
                }
            }

            template<typename Real, int BLK_X, int BLK_Y = 1, bool row_major = false>
            __global__
            void alrRcpp(
                const int ivar,
                      Real* __restrict__ out, offset_t out_stride,
                const Real* __restrict__   x, offset_t x_stride,
                int rows, int cols)
            {
                const int ivar0 = ivar - 1;

                if constexpr (row_major) {
                    const int col = blockDim.x * blockIdx.x + threadIdx.x;
                    if ((size_t)col >= cols) return;

                    for (size_t r = 0; r < rows; ++r) {
                        Real num = propr::math::log_t(x[r * x_stride + col]);
                        Real den = propr::math::log_t(x[r * x_stride + ivar0]);
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
                                Real num = propr::math::log_t(x[r + col   * x_stride]);
                                Real den = propr::math::log_t(x[r + ivar0 * x_stride]);
                                out[r + col * out_stride] = (num - den);
                            }
                        }
                        __syncthreads();
                    }
                }
            }
            
            template <typename Real, class Config>
            __global__
            //__launch_bounds__(Config::TILE * Config::BLK_N, 1, 1)
            void symRcpp(      Real* __restrict__ out, offset_t out_stride,
                         const Real* __restrict__   x, offset_t x_stride,
                         int rows, int cols){
                
                int r0 = blockIdx.x * Config::TILE + threadIdx.x;
                int c0 = blockIdx.y * Config::TILE + threadIdx.y;

                for (int dj = 0; dj < Config::TILE; dj += Config::BLK_N) {
                    int r = r0;
                    int c = c0 + dj;

                    if (r < rows && c < cols) {
                        bool can_mirror = (r < c) && (c < rows) && (r < cols);
                        int idx_rc = r + c * x_stride;
                        int idx_cr = c + r * x_stride;
                        out[r + c * out_stride] = can_mirror ? x[idx_cr] : x[idx_rc];
                    }
                }
            };

            template <typename Real, class Config>
            __global__
            //__launch_bounds__((Config::BLK_M / Config::TH_X)* (Config::BLK_M / Config::TH_Y), 1, 1)
            void phiRcpp(const bool sym,
                        Real* __restrict__ out, offset_t out_stride,
                        const Real* __restrict__   x, offset_t x_stride,
                              Real* __restrict__ row_sums,
                              Real* __restrict__ mu_sum,
                              int rows, int cols)
            {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;
                const int K = cols;
                #pragma nv_diag_default 177

                const Real* A = x;
                const Real* B = x;
                Real* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X; 
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y; 
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real S [Config::TH_Y][Config::TH_X] = {Real(0)};
                Real mu[Config::TH_Y][Config::TH_X] = {Real(0)};
                
                __syncthreads();

                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / Lanes;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / Lanes;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * Lanes;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * Lanes;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const Real* A_base = &A[(Config::BLK_M * by) * x_stride];
                const Real* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;


                auto ld_or_zero = [](const Real* __restrict__ p, int r, int c, int ld, int max_r, int max_c) {
                    return (r < max_r && c < max_c) ? propr::math::log_t(p[OFFSET(r, c, ld)]) : Real(0);
                };


                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m  = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        As[0][A_TILE_COL + lane][row_m] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                    }
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        Bs[0][row_k][col_m + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                    }
                }
                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m  = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l      = (i / A_TILE_ROW_STRIDE) * Lanes;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_a_reg[l + lane] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l      = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k  = tile_idx + B_TILE_ROW_START + i;
                            const int col_m  = B_TILE_COL;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_b_reg[l + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                            }
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem       = K - tile_base;
                    const int k_tile    = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max     = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }

                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const int k_cur = tile_base + (j + 1);

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_cur);
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                    + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l     = (i / A_TILE_ROW_STRIDE) * Lanes;
                            const int row_m = A_TILE_ROW_START + i;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][A_TILE_COL + lane][row_m] = ldg_a_reg[l + lane];
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l     = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k = B_TILE_ROW_START + i;
                            thread::store<Config::StoreModifer, Wide>(&Bs[write_stage_idx][row_k][B_TILE_COL], thread::load<Config::LoadModifer, Wide>(&ldg_b_reg[l]));
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const int k_tail  = tile_base + k_tile;
                        const int last_buf = (j_max & 1);
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_tail);
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                    + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }

                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                const Real inv_denom = (K > 1 ? Real(1) / static_cast<Real>(K - 1) : Real(1));
                PROPR_UNROLL
                for (int yy = 0; yy < Config::TH_Y; ++yy) {
                    PROPR_UNROLL
                    for (int xx = 0; xx < Config::TH_X; ++xx) {
                        S[yy][xx] *= inv_denom;
                    }
                }

                auto accum_row_mu = [&](int r, int c, Real mval){
                    if (r < M && c < M) {
                        atomicAdd(&row_sums[r], mval);
                        atomicAdd(mu_sum,       mval);
                    }
                };

                PROPR_UNROLL
                for (int local_row = 0; local_row < Config::TH_Y; ++local_row) {
                    const int row_offset = (local_row < Config::TH_Y / 2)
                        ? local_row
                        : Config::BLK_M / 2 + (local_row - Config::TH_Y / 2);
                    const int r = r0 + row_offset;

                    PROPR_UNROLL
                    for (int local_col = 0; local_col < Config::TH_X; ++local_col) {
                        const int col_offset = (local_col < Config::TH_X / 2)
                            ? local_col
                            : Config::BLK_M / 2 + (local_col - Config::TH_X / 2);
                        const int c = c0 + col_offset;
                        accum_row_mu(r, c, S[local_row][local_col]);
                    }
                }

                // --- global barrier so all reductions are complete
                auto grid = cooperative_groups::this_grid(); 
                grid.sync();

                // --- compute v_i from reductions
                const Real invM   = (M > 0 ? Real(1) / static_cast<Real>(M) : Real(0));
                const Real mu_val = (*mu_sum) * (Real(1) / (static_cast<Real>(M) * static_cast<Real>(M)));
                auto v_of = [&](int idx)->Real {
                    Real rsum = (idx < M ? row_sums[idx] : Real(0));
                    return rsum * invM - Real(0.5) * mu_val;
                };

                PROPR_UNROLL
                for (int local_row = 0; local_row < Config::TH_Y; ++local_row) {
                    const int row_offset = (local_row < Config::TH_Y / 2)
                        ? local_row
                        : Config::BLK_M / 2 + (local_row - Config::TH_Y / 2);
                    const int r = r0 + row_offset;
                    if (r >= M) continue;

                    PROPR_UNROLL
                    for (int local_col = 0; local_col < Config::TH_X; ++local_col) {
                        const int col_offset = (local_col < Config::TH_X / 2)
                            ? local_col
                            : Config::BLK_M / 2 + (local_col - Config::TH_X / 2);
                        const int c = c0 + col_offset;
                        Real val = Real(0);
                        if (c < M && r != c) {
                            const Real denom = !sym ? v_of(c) : v_of(r < c ? r : c);
                            val = S[local_row][local_col] / denom;
                        }
                        if (c < M) {
                            C[OFFSET(r, c, out_stride)] = val;
                        }
                    }
                }
            }


            template <typename Real, class Config>
            __global__
            void phiRcpp_stage1(const bool sym,
                                Real* __restrict__ out, offset_t out_stride,
                                Real* __restrict__ x,   offset_t x_stride,
                                Real* __restrict__ row_sums,
                                Real* __restrict__ mu_sum,
                                int rows, int cols)
            {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;
                const int K = cols;
                #pragma nv_diag_default 177

                const Real* A = x;
                const Real* B = x;
                Real* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X; 
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y; 
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real S [Config::TH_Y][Config::TH_X] = {Real(0)};
                Real mu[Config::TH_Y][Config::TH_X] = {Real(0)};

                __syncthreads();

                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / Lanes;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / Lanes;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * Lanes;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * Lanes;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                const Real* A_base = &A[(Config::BLK_M * by) * x_stride];
                const Real* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;

                auto ld_or_zero = [](const Real* __restrict__ p, int r, int c, int ld, int max_r, int max_c) {
                    return (r < max_r && c < max_c) ? propr::math::log_t(p[OFFSET(r, c, ld)]) : Real(0);
                };

                // --- load first A/B tiles into shared memory ---
                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m  = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        As[0][A_TILE_COL + lane][row_m] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                    }
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        Bs[0][row_k][col_m + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                    }
                }
                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base],
                        thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }
                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base],
                        thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m  = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l      = (i / A_TILE_ROW_STRIDE) * Lanes;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_a_reg[l + lane] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l      = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k  = tile_idx + B_TILE_ROW_START + i;
                            const int col_m  = B_TILE_COL;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_b_reg[l + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                            }
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem       = K - tile_base;
                    const int k_tile    = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max     = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base],
                                thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base],
                                thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const int k_cur = tile_base + (j + 1);

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_cur);
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                        + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l     = (i / A_TILE_ROW_STRIDE) * Lanes;
                            const int row_m = A_TILE_ROW_START + i;

                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][A_TILE_COL + lane][row_m] = ldg_a_reg[l + lane];
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l     = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k = B_TILE_ROW_START + i;
                            thread::store<Config::StoreModifer, Wide>(&Bs[write_stage_idx][row_k][B_TILE_COL],
                                thread::load<Config::LoadModifer, Wide>(&ldg_b_reg[l]));
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const int k_tail   = tile_base + k_tile;
                        const int last_buf = (j_max & 1);
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_tail);
                                S [thread_y][thread_x] = S[thread_y][thread_x]
                                                        + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base],
                                thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base],
                                thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                // --- finalize variance S_ij = S / (K-1) and accumulate reductions ---
                const Real inv_denom = (K > 1 ? Real(1) / static_cast<Real>(K - 1) : Real(1));
                PROPR_UNROLL
                for (int yy = 0; yy < Config::TH_Y; ++yy) {
                    PROPR_UNROLL
                    for (int xx = 0; xx < Config::TH_X; ++xx) {
                        S[yy][xx] *= inv_denom;
                    }
                }

                auto accum_row_mu = [&](int r, int c, Real mval){
                    if (r < M && c < M) {
                        atomicAdd(&row_sums[r], mval);
                        atomicAdd(mu_sum,       mval);
                    }
                };

                PROPR_UNROLL
                for (int local_row = 0; local_row < Config::TH_Y; ++local_row) {
                    const int row_offset = (local_row < Config::TH_Y / 2)
                        ? local_row
                        : Config::BLK_M / 2 + (local_row - Config::TH_Y / 2);
                    const int r = r0 + row_offset;

                    PROPR_UNROLL
                    for (int local_col = 0; local_col < Config::TH_X; ++local_col) {
                        const int col_offset = (local_col < Config::TH_X / 2)
                            ? local_col
                            : Config::BLK_M / 2 + (local_col - Config::TH_X / 2);
                        const int c = c0 + col_offset;
                        accum_row_mu(r, c, S[local_row][local_col]);
                        if (r < M && c < M) {
                            C[OFFSET(r, c, out_stride)] = S[local_row][local_col];
                        }
                    }
                }
            }


            template <typename Real, class Config>
            __global__
            void phiRcpp_stage2(const bool sym,
                                Real*  __restrict__ out, offset_t out_stride,
                                Real* __restrict__ row_sums,
                                Real* __restrict__ mu_sum,
                                int rows /* M */)
            {
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;
                const int M = rows;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X; 
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y; 
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                const int warp_id     = tid / 32;
                const int lane_id     = tid % 32;

                const int c_block_row = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int c_block_col = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                Real* C = out;

                const Real invM   = (M > 0 ? Real(1) / static_cast<Real>(M) : Real(0));
                const Real mu_val = (*mu_sum) / (static_cast<Real>(M) * static_cast<Real>(M));

                auto v_of = [&](int idx)->Real {
                    Real rsum = (idx < M ? row_sums[idx] : Real(0));
                    return rsum * invM - Real(0.5) * mu_val;
                };

                PROPR_UNROLL
                for (int local_row = 0; local_row < Config::TH_Y; ++local_row) {
                    const int row_offset = (local_row < Config::TH_Y / 2)
                        ? local_row
                        : Config::BLK_M / 2 + (local_row - Config::TH_Y / 2);
                    const int r = r0 + row_offset;
                    if (r >= M) continue;

                    PROPR_UNROLL
                    for (int local_col = 0; local_col < Config::TH_X; ++local_col) {
                        const int col_offset = (local_col < Config::TH_X / 2)
                            ? local_col
                            : Config::BLK_M / 2 + (local_col - Config::TH_X / 2);
                        const int c = c0 + col_offset;
                        Real val = Real(0);

                        if (c < M && r != c) {
                            const Real denom = !sym ? v_of(c) : v_of(r < c ? r : c);
                            if (denom != Real(0)) {
                                val = C[OFFSET(r, c, out_stride)] / denom;
                            }
                        }
                        if (c < M) {
                            C[OFFSET(r, c, out_stride)] = val;
                        }
                    }
                }
            }


            template <typename Real, class Config>
            __global__
            void 
            vlrRcpp(Real* __restrict__ out,  offset_t out_stride,
                        Real* __restrict__ x, offset_t x_stride,
                        int rows, int cols)
            {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;  // true M (nfeats)
                const int K = cols;  // true K (samples)
                #pragma nv_diag_default 177

                Real* A = x;
                Real* B = x;
                Real* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;

                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X; 
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y; 
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real S [Config::TH_Y][Config::TH_X] = {Real(0)};
                Real mu[Config::TH_Y][Config::TH_X] = {Real(0)};

                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int A_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int B_THREADS_PER_ROW = Config::BLK_K / Lanes;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_THREADS_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_THREADS_PER_ROW;

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;

                Real* A_base = A + (by * Config::BLK_M) * x_stride;
                Real* B_base = B + (bx * Config::BLK_M) * x_stride;

                {
                    const int a_m0 = tid / A_THREADS_PER_ROW;
                    const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                        const int m   = a_m0 + i;
                        const int idx = (i / A_TILE_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[idx], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            As[0][a_k + lane][m] = propr::math::log_t(ldg_a_reg[idx + lane]);
                        }
                    }
                }
                {
                    const int b_n0 = tid / B_THREADS_PER_ROW;
                    const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += B_TILE_ROW_STRIDE) {
                        const int n   = b_n0 + i;
                        const int idx = (i / B_TILE_ROW_STRIDE) * Lanes;
                        thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[idx], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k, x_stride)]));
                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            Bs[0][b_k + lane][n] = propr::math::log_t(ldg_b_reg[idx + lane]);
                        }
                    }
                }
                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }
                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;

                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int m   = a_m0 + i;
                            const int idx = (i / A_TILE_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[idx], thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, a_k + tile_idx, x_stride)]));
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_TILE_ROW_STRIDE) {
                            const int n   = b_n0 + i;
                            const int idx = (i / B_TILE_ROW_STRIDE) * Lanes;
                            thread::store<Config::StoreModifer, Wide>(&ldg_b_reg[idx], thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(n, b_k + tile_idx, x_stride)]));
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;
                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem       = K - tile_base;
                    const int k_tile    = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max     = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const int k_cur = tile_base + (j + 1);

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_cur);
                                S [thread_y][thread_x] = S[thread_y][thread_x] + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        const int a_m0 = tid / A_THREADS_PER_ROW;
                        const int a_k = (tid % A_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int m   = a_m0 + i;
                            const int idx = (i / A_TILE_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][a_k + lane][m] = propr::math::log_t(ldg_a_reg[idx + lane]);
                            }
                        }

                        const int b_n0 = tid / B_THREADS_PER_ROW;
                        const int b_k = (tid % B_THREADS_PER_ROW) * Lanes;

                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += B_TILE_ROW_STRIDE) {
                            const int n   = b_n0 + i;
                            const int idx = (i / B_TILE_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                Bs[write_stage_idx][b_k + lane][n] = propr::math::log_t(ldg_b_reg[idx + lane]);
                            }
                        }

                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const int k_tail = tile_base + k_tile;
                        const int last_buf = (j_max & 1);
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                const Real xval = a - b;

                                const Real old_mu = mu[thread_y][thread_x];
                                mu[thread_y][thread_x] = old_mu + (xval - old_mu) / static_cast<Real>(k_tail);
                                S [thread_y][thread_x] = S[thread_y][thread_x] + (xval - mu[thread_y][thread_x]) * (xval - old_mu);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const Real denom = (K > 1 ? static_cast<Real>(K - 1) : Real(1));

                PROPR_UNROLL
                for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                    const int row_offset = (thread_y < Config::TH_Y / 2)
                        ? thread_y
                        : Config::BLK_M / 2 + (thread_y - Config::TH_Y / 2);
                    const int row = by * Config::BLK_M + c_block_row + row_offset;

                    PROPR_UNROLL
                    for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                        const int col_offset = (thread_x < Config::TH_X / 2)
                            ? thread_x
                            : Config::BLK_M / 2 + (thread_x - Config::TH_X / 2);
                        const int col = bx * Config::BLK_M + c_block_col + col_offset;
                        thread::store<Config::StoreModifer, Real>(&C[OFFSET(row, col, out_stride)], S[thread_y][thread_x] / denom);
                    }
                }
            }


            template<typename Real, class Config>
            __global__
            void rhoRcpp( int ivar,
                          Real* __restrict__  out, offset_t out_stride,
                    const Real* __restrict__    x, offset_t   x_stride,   // (M_pad x K_pad)
                    const Real* __restrict__   lr, offset_t  lr_stride,   // (M_pad x K_pad)
                    int rows, int cols )                                   // rows=M_pad, cols=K(true)
            {
                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "BLK_M % TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "BLK_M % TH_X == 0");

                #pragma nv_diag_suppress 177
                const int M = rows;       // padded
                const int K = cols;       // true K
                #pragma nv_diag_default 177

                const Real* A = x;       // log(X) is applied during load
                const Real* B = x;       // same
                const Real* LR = lr;

                Real* C = out;

                const int bx = blockIdx.x, by = blockIdx.y;
                const int tx = threadIdx.x, ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK   = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK   = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As  [Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs  [Config::BLK_K][Config::BLK_M];
                __shared__ Real LR_A[Config::BLK_K][Config::BLK_M];
                __shared__ Real LR_B[Config::BLK_K][Config::BLK_M];

                Real S [Config::TH_Y][Config::TH_X] = {Real(0)};
                Real mu[Config::TH_Y][Config::TH_X] = {Real(0)};
                Real mu_lr_i[Config::TH_Y] = {Real(0)}, mu_lr_j[Config::TH_X] = {Real(0)};
                Real S_lr_i [Config::TH_Y] = {Real(0)}, S_lr_j [Config::TH_X] = {Real(0)};

                Real frag_xi [2][Config::TH_Y];
                Real frag_xj [2][Config::TH_X];
                Real frag_lr_i[2][Config::TH_Y];
                Real frag_lr_j[2][Config::TH_X];

                const int A_THREADS_PER_ROW = Config::BLK_K / Lanes;
                const int ROW_STRIDE        = THREAD_NUM_PER_BLOCK / A_THREADS_PER_ROW;

                const Real* A_base    = &A[(Config::BLK_M * by) * x_stride];
                const Real* B_base    = &B[(Config::BLK_M * bx) * x_stride];
                const Real* LR_A_base = &LR[(Config::BLK_M * by) * lr_stride];
                const Real* LR_B_base = &LR[(Config::BLK_M * bx) * lr_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;

                PROPR_NO_UNROLL
                for (int tile_base = 0; tile_base < K; tile_base += Config::BLK_K) {
                    const int k_tile = min(Config::BLK_K, K - tile_base);

                    const int m0  = tid / A_THREADS_PER_ROW;
                    const int k_base  = (tid % A_THREADS_PER_ROW) * Lanes;

                    PROPR_UNROLL
                    for (int i = 0; i < Config::BLK_M; i += ROW_STRIDE) {
                        const int m = m0 + i;
                        Wide xa = thread::load<Config::LoadModifer, Wide>(&A_base[OFFSET(m, tile_base + k_base, x_stride)]);
                        Wide xb = thread::load<Config::LoadModifer, Wide>(&B_base[OFFSET(m, tile_base + k_base, x_stride)]);
                        Wide la = thread::load<Config::LoadModifer, Wide>(&LR_A_base[OFFSET(m, tile_base + k_base, lr_stride)]);
                        Wide lb = thread::load<Config::LoadModifer, Wide>(&LR_B_base[OFFSET(m, tile_base + k_base, lr_stride)]);

                        PROPR_UNROLL
                        for (int lane = 0; lane < Lanes; ++lane) {
                            As  [k_base + lane][m] = propr::math::log_t(propr::lane_at(xa, lane));
                            Bs  [k_base + lane][m] = propr::math::log_t(propr::lane_at(xb, lane));
                            LR_A[k_base + lane][m] = propr::lane_at(la, lane);
                            LR_B[k_base + lane][m] = propr::lane_at(lb, lane);
                        }
                    }

                    __syncthreads();

                    if (k_tile > 0) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_xi[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][a_tile_index + offset]));
                            thread::store<Config::StoreModifer, Wide>(&frag_lr_i[0][base], thread::load<Config::LoadModifer, Wide>(&LR_A[0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_xj[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][b_tile_index + offset]));
                            thread::store<Config::StoreModifer, Wide>(&frag_lr_j[0][base], thread::load<Config::LoadModifer, Wide>(&LR_B[0][b_tile_index + offset]));
                        }
                    }

                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_NO_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_xi[(j + 1)&1][base], thread::load<Config::LoadModifer, Wide>(&As[j + 1][a_tile_index + offset]));
                            thread::store<Config::StoreModifer, Wide>(&frag_lr_i[(j + 1)&1][base], thread::load<Config::LoadModifer, Wide>(&LR_A[j + 1][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_xj[(j + 1)&1][base], thread::load<Config::LoadModifer, Wide>(&Bs[j + 1][b_tile_index + offset]));
                            thread::store<Config::StoreModifer, Wide>(&frag_lr_j[(j + 1)&1][base], thread::load<Config::LoadModifer, Wide>(&LR_B[j + 1][b_tile_index + offset]));
                        }

                        const int k_cur = tile_base + (j + 1);

                        PROPR_UNROLL
                        for (int ry = 0; ry < Config::TH_Y; ++ry) {
                            const Real a = frag_xi[j & 1][ry];
                            PROPR_UNROLL
                            for (int rx = 0; rx < Config::TH_X; ++rx) {
                                const Real b    = frag_xj[j & 1][rx];
                                const Real xval = a - b;
                                const Real old_mu = mu[ry][rx];
                                mu[ry][rx] = old_mu + (xval - old_mu) / static_cast<Real>(k_cur);
                                S [ry][rx] = S[ry][rx] + (xval - mu[ry][rx]) * (xval - old_mu);
                            }
                        }

                        PROPR_UNROLL
                        for (int ry = 0; ry < Config::TH_Y; ++ry) {
                            const Real vi_old = mu_lr_i[ry];
                            const Real vi_val = frag_lr_i[j & 1][ry];
                            mu_lr_i[ry] = vi_old + (vi_val - vi_old) / static_cast<Real>(k_cur);
                            S_lr_i [ry] = S_lr_i[ry] + (vi_val - mu_lr_i[ry]) * (vi_val - vi_old);
                        }
                        PROPR_UNROLL
                        for (int rx = 0; rx < Config::TH_X; ++rx) {
                            const Real vj_old = mu_lr_j[rx];
                            const Real vj_val = frag_lr_j[j & 1][rx];
                            mu_lr_j[rx] = vj_old + (vj_val - vj_old) / static_cast<Real>(k_cur);
                            S_lr_j [rx] = S_lr_j[rx] + (vj_val - mu_lr_j[rx]) * (vj_val - vj_old);
                        }
                    }

                    if (k_tile > 0) {
                        const int k_tail  = tile_base + k_tile;
                        const int last_buf = (j_max & 1);
                        PROPR_UNROLL
                        for (int ry = 0; ry < Config::TH_Y; ++ry) {
                            const Real a = frag_xi[last_buf][ry];
                            PROPR_UNROLL
                            for (int rx = 0; rx < Config::TH_X; ++rx) {
                                const Real b    = frag_xj[last_buf][rx];
                                const Real xval = a - b;
                                const Real old_mu = mu[ry][rx];
                                mu[ry][rx] = old_mu + (xval - old_mu) / static_cast<Real>(k_tail);
                                S [ry][rx] = S[ry][rx] + (xval - mu[ry][rx]) * (xval - old_mu);
                            }
                        }
                        for (int ry = 0; ry < Config::TH_Y; ++ry) {
                            const Real vi_old = mu_lr_i[ry];
                            const Real vi_val = frag_lr_i[last_buf][ry];
                            mu_lr_i[ry] = vi_old + (vi_val - vi_old) / static_cast<Real>(k_tail);
                            S_lr_i [ry] = S_lr_i[ry] + (vi_val - mu_lr_i[ry]) * (vi_val - vi_old);
                        }
                        PROPR_UNROLL
                        for (int rx = 0; rx < Config::TH_X; ++rx) {
                            const Real vj_old = mu_lr_j[rx];
                            const Real vj_val = frag_lr_j[last_buf][rx];
                            mu_lr_j[rx] = vj_old + (vj_val - vj_old) / static_cast<Real>(k_tail);
                            S_lr_j [rx] = S_lr_j[rx] + (vj_val - mu_lr_j[rx]) * (vj_val - vj_old);
                        }
                    }

                    __syncthreads();
                }

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const Real denom = (K > 1 ? static_cast<Real>(K - 1) : Real(1));

                const int r0 = Config::BLK_M * by + c_block_row;
                const int c0 = Config::BLK_M * bx + c_block_col;

                Real var_i_rows[Config::TH_Y];
                Real var_j_cols[Config::TH_X];
                PROPR_UNROLL
                for (int i = 0; i < Config::TH_Y; ++i) {
                    var_i_rows[i] = S_lr_i[i] / denom;
                }
                PROPR_UNROLL
                for (int j = 0; j < Config::TH_X; ++j) {
                    var_j_cols[j] = S_lr_j[j] / denom;
                }

                const int ivar0 = ivar - 1;

                auto rho_value = [&](int local_row, int local_col, int global_row, int global_col) -> Real {
                    const Real is_ivar = static_cast<Real>((global_row == ivar0) | (global_col == ivar0));
                    const Real is_diag_ivar = static_cast<Real>((global_row == ivar0) & (global_col == ivar0));
                    return (Real(1) - ((S[local_row][local_col] / denom) / (var_i_rows[local_row] + var_j_cols[local_col])))
                        * (Real(1) - is_ivar) + is_diag_ivar;
                };

                PROPR_UNROLL
                for (int local_row = 0; local_row < Config::TH_Y; ++local_row) {
                    const int row_offset = (local_row < Config::TH_Y / 2)
                        ? local_row
                        : Config::BLK_M / 2 + (local_row - Config::TH_Y / 2);
                    const int global_row = r0 + row_offset;

                    PROPR_UNROLL
                    for (int local_col = 0; local_col < Config::TH_X; ++local_col) {
                        const int col_offset = (local_col < Config::TH_X / 2)
                            ? local_col
                            : Config::BLK_M / 2 + (local_col - Config::TH_X / 2);
                        const int global_col = c0 + col_offset;
                        thread::store<Config::StoreModifer, Real>(
                            &C[OFFSET(global_row, global_col, out_stride)],
                            rho_value(local_row, local_col, global_row, global_col)
                        );
                    }
                }
            }

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


            template <typename Real, class Config>
            __global__
            void 
            linRcpp(Real* __restrict__ out, offset_t out_stride,
                    Real* __restrict__ rho, int rho_stride,
                    Real* __restrict__ x, offset_t x_stride,
                    int rows, int cols){

                using Wide = propr::cuda_wide_vector_t<Real>;
                constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

                static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
                static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
                static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

                const int M = rows;      // features (nfeats)
                const int K = cols;      // samples  (N_samples)

                Real* A = x;
                Real* B = x;
                Real* C = out;

                const int bx = blockIdx.x;
                const int by = blockIdx.y;
                const int tx = threadIdx.x;
                const int ty = threadIdx.y;

                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;
                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
                __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

                Real Sa[Config::TH_Y] = {Real(0)};
                Real Sb[Config::TH_X] = {Real(0)};
                Real mu_a[Config::TH_Y] = {Real(0)};
                Real mu_b[Config::TH_X] = {Real(0)};
                Real accum[Config::TH_Y][Config::TH_X] = {Real(0)};

                Real frag_a[2][Config::TH_Y];
                Real frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * Lanes);
                Real ldg_a_reg[Lanes * ldg_num_a];
                Real ldg_b_reg[Lanes * ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / Lanes;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / Lanes;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = (tid % A_TILE_THREAD_PER_ROW) * Lanes;
                const int B_TILE_COL = (tid % B_TILE_THREAD_PER_ROW) * Lanes;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                Real* A_base = &A[(Config::BLK_M * by) * x_stride];
                Real* B_base = &B[(Config::BLK_M * bx) * x_stride];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index = (warp_id / 2) * 16 + (lane_id / 8) * 4;
                const int b_tile_index = (warp_id % 2) * 32 + (lane_id % 8) * 4;

                auto ld_or_zero = [](Real* __restrict__ p,
                                    int r, int c, int ld, int max_r, int max_c) {
                    return (r >= 0 && r < max_r && c >= 0 && c < max_c) ? p[OFFSET(r, c, ld)] : Real(0);
                };

                auto rho_at = [&](int r, int c) -> Real {
                    return ld_or_zero(rho, r, c, rho_stride, M, M);
                };

                auto store_if_in_bounds = [&](int r, int c, Real v) {
                    if (r < M && c < M) C[OFFSET(r, c, out_stride)] = v;
                };

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                    const int row_m = A_TILE_ROW_START + i;
                    const int base_k = A_TILE_COL;

                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        As[0][A_TILE_COL + lane][row_m] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                    }
                }

                PROPR_UNROLL
                for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    const int row_k = B_TILE_ROW_START + i;
                    const int col_m = B_TILE_COL;
                    PROPR_UNROLL
                    for (int lane = 0; lane < Lanes; ++lane) {
                        Bs[0][row_k][col_m + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                    }
                }
                __syncthreads();

                PROPR_UNROLL
                for (int base = 0; base < Config::TH_Y; base += Lanes) {
                    const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[0][0][a_tile_index + offset]));
                }
                PROPR_UNROLL
                for (int base = 0; base < Config::TH_X; base += Lanes) {
                    const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                    thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[0][0][b_tile_index + offset]));
                }

                int write_stage_idx = 1;
                int tile_idx = 0;

                do {
                    tile_idx += Config::BLK_K;
                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int row_m = A_TILE_ROW_START + i;
                            const int base_k = A_TILE_COL + tile_idx;
                            const int l = (i / A_TILE_ROW_STRIDE) * Lanes;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_a_reg[l + lane] = ld_or_zero(A_base, row_m, base_k + lane, x_stride, M, K);
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k = tile_idx + B_TILE_ROW_START + i;
                            const int col_m = B_TILE_COL;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                ldg_b_reg[l + lane] = ld_or_zero(B_base, col_m + lane, row_k, x_stride, M, K);
                            }
                        }
                    }

                    const int load_stage_idx = write_stage_idx ^ 1;

                    const int tile_base = tile_idx - Config::BLK_K;
                    const int rem = K - tile_base;
                    const int k_tile = (rem < Config::BLK_K ? rem : Config::BLK_K);
                    const int j_max = (k_tile > 0 ? k_tile - 1 : 0);

                    PROPR_UNROLL
                    for (int j = 0; j < j_max; ++j) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j + 1)][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[(j + 1) & 1][base], thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j + 1)][b_tile_index + offset]));
                        }

                        const Real n = static_cast<Real>(tile_base + (j + 1));

                        // Update mu_b, Sb
                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[j & 1][thread_x];
                            Real db = b - mu_b[thread_x];
                            Real mu_b_new = mu_b[thread_x] + db / n;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        // Update mu_a, Sa, and cross
                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[j & 1][thread_y];
                            Real da = a - mu_a[thread_y];
                            Real mu_a_new = mu_a[thread_y] + da / n;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[j & 1][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_M; i += A_TILE_ROW_STRIDE) {
                            const int l = (i / A_TILE_ROW_STRIDE) * Lanes;
                            const int row_m = A_TILE_ROW_START + i;
                            PROPR_UNROLL
                            for (int lane = 0; lane < Lanes; ++lane) {
                                As[write_stage_idx][A_TILE_COL + lane][row_m] = ldg_a_reg[l + lane];
                            }
                        }
                        PROPR_UNROLL
                        for (int i = 0; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            const int l = (i / B_TILE_ROW_STRIDE) * Lanes;
                            const int row_k = B_TILE_ROW_START + i;
                            thread::store<Config::StoreModifer, Wide>(&Bs[write_stage_idx][row_k][B_TILE_COL], thread::load<Config::LoadModifer, Wide>(&ldg_b_reg[l]));
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    if (k_tile > 0) {
                        const Real n_tail = static_cast<Real>(tile_base + k_tile);
                        const int last_buf = (j_max & 1);

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            const Real b = frag_b[last_buf][thread_x];
                            Real db = b - mu_b[thread_x];
                            Real mu_b_new = mu_b[thread_x] + db / n_tail;
                            Sb[thread_x] += db * (b - mu_b_new);
                            mu_b[thread_x] = mu_b_new;
                        }

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            const Real a = frag_a[last_buf][thread_y];
                            Real da = a - mu_a[thread_y];
                            Real mu_a_new = mu_a[thread_y] + da / n_tail;
                            Sa[thread_y] += da * (a - mu_a_new);
                            mu_a[thread_y] = mu_a_new;

                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                const Real b = frag_b[last_buf][thread_x];
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            }
                        }
                    }

                    if (tile_idx < K) {
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_Y; base += Lanes) {
                            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base], thread::load<Config::LoadModifer, Wide>(&As[(write_stage_idx ^ 1)][0][a_tile_index + offset]));
                        }
                        PROPR_UNROLL
                        for (int base = 0; base < Config::TH_X; base += Lanes) {
                            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base], thread::load<Config::LoadModifer, Wide>(&Bs[(write_stage_idx ^ 1)][0][b_tile_index + offset]));
                        }
                    }
                } while (tile_idx < K);

                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                const Real eps     = Real(1e-20);
                const Real eps_rho = Real(1e-7); // clamp for atanh

                auto atanh_clamped = [&](Real x) -> Real {
                    Real xc = propr::math::clamp_t(x, Real(-1) + eps_rho, Real(1) - eps_rho);
                    return Real(0.5) * propr::math::log_t((Real(1) + xc) / (Real(1) - xc));
                };


                Real Creg[Config::TH_Y][Config::TH_X];

                PROPR_UNROLL
                for (int i = 0; i < Config::TH_Y; ++i) {
                    const int row_offset = (i < Config::TH_Y / 2)
                        ? i
                        : Config::BLK_M / 2 + (i - Config::TH_Y / 2);
                    const int rr = Config::BLK_M * by + c_block_row + row_offset;
                    const Real inva = (Sa[i] > eps) ? propr::math::rsqrt_t(Sa[i]) : Real(0);

                    PROPR_UNROLL
                    for (int j = 0; j < Config::TH_X; ++j) {
                        const int col_offset = (j < Config::TH_X / 2)
                            ? j
                            : Config::BLK_M / 2 + (j - Config::TH_X / 2);
                        const int cc = Config::BLK_M * bx + c_block_col + col_offset;

                        const Real invb = (Sb[j] > eps) ? propr::math::rsqrt_t(Sb[j]) : Real(0);
                        Real r = accum[i][j] * inva * invb;
                        r = propr::math::clamp_t(r, Real(-1), Real(1));

                        Real out_val = Real(1);
                        const Real rho_ij = rho_at(rr, cc);
                        if (rr > cc) {
                            // lower triangle: z = atanh(rho)
                            out_val = atanh_clamped(rho_ij);
                        } else  if (rr < cc) {
                            // upper triangle: variance formula
                            const Real r2     = propr::math::max_t(r * r, eps);
                            const Real rho2   = rho_ij * rho_ij;
                            const Real denom  = propr::math::max_t((Real(1) - rho2) * r2 * (static_cast<Real>(K) - Real(2)), eps);
                            const Real num    = (Real(1) - r2) * rho2;
                            out_val = num / denom;
                        }

                        Creg[i][j] = out_val;
                    }
                }

                PROPR_UNROLL
                for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                    const int row_offset = (thread_y < Config::TH_Y / 2)
                        ? thread_y
                        : Config::BLK_M / 2 + (thread_y - Config::TH_Y / 2);
                    const int row = Config::BLK_M * by + c_block_row + row_offset;

                    PROPR_UNROLL
                    for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                        const int col_offset = (thread_x < Config::TH_X / 2)
                            ? thread_x
                            : Config::BLK_M / 2 + (thread_x - Config::TH_X / 2);
                        const int col = Config::BLK_M * bx + c_block_col + col_offset;
                        store_if_in_bounds(row, col, Creg[thread_y][thread_x]);
                    }
                }
            }


            template <typename Real>
            __global__
            void lltRcpp(
                    Real * out,
                    size_t nfeats,
                    Real * __restrict__ X,
                    offset_t x_stride)
            {
                const offset_t total_pairs = (nfeats * (nfeats - 1)) / 2;

                for (offset_t k = blockDim.x * blockIdx.x + threadIdx.x;
                    k < total_pairs;
                    k += gridDim.x * blockDim.x) {

                    int i = 0;
                    int j = 0;
                    thread::lower_triangle_pair_from_index(static_cast<std::size_t>(k), i, j);
                    out[k] = X[j * x_stride + i];
                }
            }


            template <typename Real>
            __global__
            void urtRcpp(
                     Real * __restrict__ out, size_t n,
                     Real * __restrict__   X,  offset_t x_stride){
                // TODO: Check f4 perf
                offset_t total_pairs = (n * (n - 1)) / 2;

                PROPR_UNROLL
                for (offset_t k = blockDim.x * blockIdx.x + threadIdx.x; k < total_pairs; k += gridDim.x * blockDim.x) {
                    int i = 0;
                    int j = 0;
                    thread::lower_triangle_pair_from_index(static_cast<std::size_t>(k), i, j);
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
                    int i = 0;
                    int j = 0;
                    thread::lower_triangle_pair_from_index(static_cast<std::size_t>(k), i, j);

                    partner[k] = i + 1;
                    pair[k]    = j + 1;
                }
            };

            template <typename Real>
            __global__
            void half2mat(Real*       __restrict__ out,
                          offset_t out_stride,
                          const Real* __restrict__ X,
                          size_t n) {
                // TODO: Check f4 perf
                offset_t total_pairs = n * (n - 1) / 2;

                PROPR_UNROLL
                for (offset_t k = blockIdx.x * blockDim.x + threadIdx.x; k < total_pairs; k += blockDim.x * gridDim.x) {
                    int i = 0;
                    int j = 0;
                    thread::lower_triangle_pair_from_index(static_cast<std::size_t>(k), i, j);
                    out[j * out_stride + i] = X[k];
                    out[i * out_stride + j] = X[k];
                }
            };

            template <typename Real>
            __global__
            void vector2mat(
                Real*      __restrict__ out,  offset_t    out_stride,
                const Real* __restrict__ X,
                const int* __restrict__ i_vec, const int* __restrict__ j_vec,
                size_t ni
            ){
                // TODO: Check f4 perf
                PROPR_UNROLL
                for (offset_t idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ni; idx += blockDim.x * gridDim.x) {
                    out[i_vec[idx] - 1 + (j_vec[idx]-1) * out_stride ] = X[idx];
                    out[j_vec[idx] - 1 + (i_vec[idx]-1) * out_stride ] = X[idx];
                }
            };

            template <typename Real>
            __global__
            void ratiosRcpp(
                Real*      __restrict__ out,
                offset_t    out_stride,
                const Real* __restrict__ X,
                offset_t    X_stride,
                size_t nfeats, size_t nsamps){

                const offset_t total_pairs = nfeats * (nfeats - 1) / 2;
                const offset_t total_elems = total_pairs * nsamps;

                PROPR_UNROLL
                for (offset_t idx = blockIdx.x * blockDim.x + threadIdx.x; idx < total_elems; idx += blockDim.x * gridDim.x) {
                    offset_t k = idx / nsamps;
                    offset_t s = static_cast<offset_t>(idx % nsamps);
                    int i = 0;
                    int j = 0;
                    thread::lower_triangle_pair_from_index(static_cast<std::size_t>(k), i, j);
                    out[s + k * out_stride] = X[s + i * X_stride] / X[s + j * X_stride];
                }
            };

            template <typename Real>
            __global__
            void results2matRcpp(
                      Real* __restrict__     out, offset_t     out_stride,
                const Real* __restrict__ results, offset_t results_stride,
                Real diagonal, size_t n){
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
