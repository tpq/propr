#pragma once

#include <cuda_runtime.h>

#include <cub/cub.cuh>

#include <propr/data/types.h>
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

            template<int BLK_X>
            __global__
            void centerNumericMatrix(
                     float * out,
                     offset_t out_stride,
                     float * __restrict__ x, 
                     offset_t x_stride,
                     size_t rows, size_t cols) {
                // TODO: impl a col major variant
                // TODO: investigate pipline sol (with split warps prod-cons)
                // TODO: f4
                const int col = blockDim.x * blockIdx.x + threadIdx.x;
                if ((size_t)col >= cols) return;

                float mean = 0.0;
                for (size_t r = 0; r < rows; ++r) {
                    float v = x[r * x_stride + col];
                    mean += (v - mean) / (r + 1);
                }

                for (size_t r = 0; r < rows; ++r) {
                    float v = x[r * x_stride + col];
                    out[r * out_stride + col] = (v - mean);
                }
            };


            struct cor_config { // TODO: move into trait like system
                const static int BLK_M = 128;
                const static int BLK_K = 8;
                const static int TH_Y  = 8;
                const static int TH_X  = 8;
            };

            template <class Config>
            __global__
            void corRcpp(float*__restrict__ out, offset_t out_stride,
                         float* __restrict__ x,  offset_t x_stride,
                         int rows, int cols) {
                // this version is the on the fly algo
                //(cov a * b)/(vr_a * vr_b) 
                const int M = rows;
                const int K = cols;

                float* A  = x;
                float* B  = x;
                float* C  = x;

                int bx = blockIdx.x;
                int by = blockIdx.y;
                // if (bx > by) return;

                int tx = threadIdx.x;
                int ty = threadIdx.y;
                
                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;

                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];


                float Sa    [Config::TH_Y] = {0.0f};
                float Sb    [Config::TH_X] = {0.0f};

                float mu_a  [Config::TH_Y] = {0.0f};
                float mu_b  [Config::TH_X] = {0.0f};
                
                float accum [Config::TH_Y][Config::TH_X] = {0.0f};

                int na = 1, nb = 1;

                float frag_a[2][Config::TH_Y];
                float frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);
                float ldg_a_reg[4*ldg_num_a];
                float ldg_b_reg[4*ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = tid % A_TILE_THREAD_PER_ROW * 4;
                const int B_TILE_COL = tid % B_TILE_THREAD_PER_ROW * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                A = &A[(Config::BLK_M * by) * K];
                B = &B[(Config::BLK_M * bx) * K];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  warp_id/2*16 + lane_id/8*4;
                const int b_tile_index =  warp_id%2*32 + lane_id%8*4;
                
                PROPR_UNROLL
                for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                    int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                    FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(
                        A_TILE_ROW_START + i, 
                        A_TILE_COL, 
                        K )]);
                    As[0][A_TILE_COL  ][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index  ];
                    As[0][A_TILE_COL+1][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+1];
                    As[0][A_TILE_COL+2][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+2];
                    As[0][A_TILE_COL+3][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+3];
                }

                PROPR_UNROLL
                for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    int row_k = B_TILE_ROW_START + i;
                    int col_n = B_TILE_COL;
                    Bs[0][row_k][col_n + 0] = B[OFFSET(col_n + 0, row_k, K)];
                    Bs[0][row_k][col_n + 1] = B[OFFSET(col_n + 1, row_k, K)];
                    Bs[0][row_k][col_n + 2] = B[OFFSET(col_n + 2, row_k, K)];
                    Bs[0][row_k][col_n + 3] = B[OFFSET(col_n + 3, row_k, K)];
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

                    if(tile_idx < K){
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                            int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                            FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(A_TILE_ROW_START + i,  
                                                                                       A_TILE_COL + tile_idx,
                                                                                       K )]);
                        }
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            int ldg_index = i / B_TILE_ROW_STRIDE * 4;
                            int row_k = tile_idx + B_TILE_ROW_START + i;
                            int col_n = B_TILE_COL;

                            ldg_b_reg[ldg_index + 0] = B[OFFSET(col_n + 0, row_k, K)];
                            ldg_b_reg[ldg_index + 1] = B[OFFSET(col_n + 1, row_k, K)];
                            ldg_b_reg[ldg_index + 2] = B[OFFSET(col_n + 2, row_k, K)];
                            ldg_b_reg[ldg_index + 3] = B[OFFSET(col_n + 3, row_k, K)];
                        }
                    }

                    int load_stage_idx = write_stage_idx ^ 1;


                    PROPR_UNROLL
                    for(int j=0; j < Config::BLK_K - 1; ++j){
                        FETCH_FLOAT4(frag_a[(j+1)%2][0]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j+1)%2][4]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index + 64]);

                        FETCH_FLOAT4(frag_b[(j+1)%2][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j+1)%2][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index + 64]);
                        

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            auto a  = frag_a[j%2][thread_y];
                            auto da = a  - mu_a[thread_y];
                            mu_a[thread_y] += da / na;
                            Sa[thread_y]   += da * (a - mu_a[thread_y]); 
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                auto b = frag_b[j%2][thread_x];
                                auto db = b  - mu_b[thread_x];
                                mu_b[thread_x] += db / nb;
                                Sb  [thread_x] += db * (b - mu_b[thread_x]); 

                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                                nb+=1;
                            }
                            na +=1;
                        }
                    }

                    if(tile_idx < K) {
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                            int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                            As[write_stage_idx][A_TILE_COL  ][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index];
                            As[write_stage_idx][A_TILE_COL+1][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+1];
                            As[write_stage_idx][A_TILE_COL+2][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+2];
                            As[write_stage_idx][A_TILE_COL+3][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+3];
                        }
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            int ldg_index = i / B_TILE_ROW_STRIDE * 4;
                            FETCH_FLOAT4(Bs[write_stage_idx][B_TILE_ROW_START + i][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[ldg_index]);
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[load_stage_idx^1][0][a_tile_index]);
                    FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[load_stage_idx^1][0][a_tile_index + 64]);

                    FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[load_stage_idx^1][0][b_tile_index]);
                    FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[load_stage_idx^1][0][b_tile_index + 64]);

                    PROPR_UNROLL
                    for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {

                        auto a = frag_a[1][thread_y];
                        auto da = a  - mu_a[thread_y];
                        mu_a[thread_y] += da / na;
                        Sa[thread_y]   += da * (a - mu_a[thread_y]); 

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            auto b = frag_b[1][thread_x];
                            auto db = b  - mu_b[thread_x];
                            mu_b[thread_x] += db / nb;
                            Sb[thread_x]   += db * (b - mu_b[thread_x]); 
                            accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            nb+=1;
                        }
                        na +=1;
                    }
                } while(tile_idx < K);
                
                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;

                for (int i = 0; i < 4; ++i) {
                    float4 tmp0 = make_float4(
                        accum[i][0] / sqrt(Sa[i] * Sb[0]),
                        accum[i][1] / sqrt(Sa[i] * Sb[1]), 
                        accum[i][2] / sqrt(Sa[i] * Sb[2]), 
                        accum[i][3] / sqrt(Sa[i] * Sb[3])
                    );

                    float4 tmp1 = make_float4(
                        accum[i][4] / sqrt(Sa[i] * Sb[4]),
                        accum[i][5] / sqrt(Sa[i] * Sb[5]),
                        accum[i][6] / sqrt(Sa[i] * Sb[6]),
                        accum[i][7] / sqrt(Sa[i] * Sb[7])
                    );

                    float4 tmp2 = make_float4(
                        accum[i+4][0] / sqrt(Sa[i+4] * Sb[0]),
                        accum[i+4][1] / sqrt(Sa[i+4] * Sb[1]),
                        accum[i+4][2] / sqrt(Sa[i+4] * Sb[2]),
                        accum[i+4][3] / sqrt(Sa[i+4] * Sb[3])
                    );

                    float4 tmp3 = make_float4(
                        accum[i+4][4] / sqrt(Sa[i+4] * Sb[4]) ,
                        accum[i+4][5] / sqrt(Sa[i+4] * Sb[5]) ,
                        accum[i+4][6] / sqrt(Sa[i+4] * Sb[6]) ,
                        accum[i+4][7] / sqrt(Sa[i+4] * Sb[7]) 
                    );

                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col,      M)]) = tmp0;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col + 64, M)]) = tmp1;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col,      M)]) = tmp2;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col + 64, M)]) = tmp3;
                }
            };

            struct cov_config { // TODO: move into trait like system
                const static int BLK_M = 128;
                const static int BLK_K = 8;
                const static int TH_Y  = 8;
                const static int TH_X  = 8;
            };

            template <class Config>
            __global__
            void 
            covRcpp(
                const int norm_type,
                float*__restrict__ out, offset_t out_stride,
                float* __restrict__ x,  offset_t x_stride,
                int rows, int cols
            ){

                const int M = rows;
                const int K = cols;

                float* A  = x;
                float* B  = x;
                float* C  = x;

                int bx = blockIdx.x;
                int by = blockIdx.y;
                // if (bx > by) return;

                int tx = threadIdx.x;
                int ty = threadIdx.y;
                
                const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
                const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
                const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;

                const int tid = ty * THREAD_X_PER_BLOCK + tx;

                __shared__ float As[2][Config::BLK_K][Config::BLK_M];
                __shared__ float Bs[2][Config::BLK_K][Config::BLK_M];


                float Sa    [Config::TH_Y] = {0.0f};
                float Sb    [Config::TH_X] = {0.0f};

                float mu_a  [Config::TH_Y] = {0.0f};
                float mu_b  [Config::TH_X] = {0.0f};
                
                float accum [Config::TH_Y][Config::TH_X] = {0.0f};

                int na = 1, nb = 1;

                float frag_a[2][Config::TH_Y];
                float frag_b[2][Config::TH_X];

                const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * 4);
                const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * 4);
                float ldg_a_reg[4*ldg_num_a];
                float ldg_b_reg[4*ldg_num_b];

                const int A_TILE_THREAD_PER_ROW = Config::BLK_K / 4;
                const int B_TILE_THREAD_PER_ROW = Config::BLK_M / 4;

                const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

                const int A_TILE_COL = tid % A_TILE_THREAD_PER_ROW * 4;
                const int B_TILE_COL = tid % B_TILE_THREAD_PER_ROW * 4;

                const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
                const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

                A = &A[(Config::BLK_M * by) * K];
                B = &B[(Config::BLK_M * bx) * K];

                const int warp_id = tid / 32;
                const int lane_id = tid % 32;
                const int a_tile_index =  warp_id/2*16 + lane_id/8*4;
                const int b_tile_index =  warp_id%2*32 + lane_id%8*4;
                
                PROPR_UNROLL
                for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                    int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                    FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(
                        A_TILE_ROW_START + i, 
                        A_TILE_COL, 
                        K )]);
                    As[0][A_TILE_COL  ][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index  ];
                    As[0][A_TILE_COL+1][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+1];
                    As[0][A_TILE_COL+2][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+2];
                    As[0][A_TILE_COL+3][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+3];
                }

                PROPR_UNROLL
                for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                    int row_k = B_TILE_ROW_START + i;
                    int col_n = B_TILE_COL;
                    Bs[0][row_k][col_n + 0] = B[OFFSET(col_n + 0, row_k, K)];
                    Bs[0][row_k][col_n + 1] = B[OFFSET(col_n + 1, row_k, K)];
                    Bs[0][row_k][col_n + 2] = B[OFFSET(col_n + 2, row_k, K)];
                    Bs[0][row_k][col_n + 3] = B[OFFSET(col_n + 3, row_k, K)];
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

                    if(tile_idx < K){
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                            int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                            FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(A_TILE_ROW_START + i,  
                                                                                       A_TILE_COL + tile_idx,
                                                                                       K )]);
                        }
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            int ldg_index = i / B_TILE_ROW_STRIDE * 4;
                            int row_k = tile_idx + B_TILE_ROW_START + i;
                            int col_n = B_TILE_COL;

                            ldg_b_reg[ldg_index + 0] = B[OFFSET(col_n + 0, row_k, K)];
                            ldg_b_reg[ldg_index + 1] = B[OFFSET(col_n + 1, row_k, K)];
                            ldg_b_reg[ldg_index + 2] = B[OFFSET(col_n + 2, row_k, K)];
                            ldg_b_reg[ldg_index + 3] = B[OFFSET(col_n + 3, row_k, K)];
                        }
                    }

                    int load_stage_idx = write_stage_idx ^ 1;


                    PROPR_UNROLL
                    for(int j=0; j < Config::BLK_K - 1; ++j){
                        FETCH_FLOAT4(frag_a[(j+1)%2][0]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index]);
                        FETCH_FLOAT4(frag_a[(j+1)%2][4]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index + 64]);

                        FETCH_FLOAT4(frag_b[(j+1)%2][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index]);
                        FETCH_FLOAT4(frag_b[(j+1)%2][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index + 64]);
                        

                        PROPR_UNROLL
                        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                            auto a  = frag_a[j%2][thread_y];
                            auto da = a  - mu_a[thread_y];
                            mu_a[thread_y] += da / na;
                            PROPR_UNROLL
                            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                                auto b = frag_b[j%2][thread_x];
                                auto db = b  - mu_b[thread_x];
                                mu_b[thread_x] += db / nb;
                                accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                                nb+=1;
                            }
                            na +=1;
                        }
                    }

                    if(tile_idx < K) {
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                            int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                            As[write_stage_idx][A_TILE_COL  ][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index];
                            As[write_stage_idx][A_TILE_COL+1][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+1];
                            As[write_stage_idx][A_TILE_COL+2][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+2];
                            As[write_stage_idx][A_TILE_COL+3][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index+3];
                        }
                        PROPR_UNROLL
                        for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                            int ldg_index = i / B_TILE_ROW_STRIDE * 4;
                            FETCH_FLOAT4(Bs[write_stage_idx][B_TILE_ROW_START + i][B_TILE_COL]) = FETCH_FLOAT4(ldg_b_reg[ldg_index]);
                        }
                        __syncthreads();
                        write_stage_idx ^= 1;
                    }

                    FETCH_FLOAT4(frag_a[0][0]) = FETCH_FLOAT4(As[load_stage_idx^1][0][a_tile_index]);
                    FETCH_FLOAT4(frag_a[0][4]) = FETCH_FLOAT4(As[load_stage_idx^1][0][a_tile_index + 64]);

                    FETCH_FLOAT4(frag_b[0][0]) = FETCH_FLOAT4(Bs[load_stage_idx^1][0][b_tile_index]);
                    FETCH_FLOAT4(frag_b[0][4]) = FETCH_FLOAT4(Bs[load_stage_idx^1][0][b_tile_index + 64]);

                    PROPR_UNROLL
                    for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {

                        auto a = frag_a[1][thread_y];
                        auto da = a  - mu_a[thread_y];
                        mu_a[thread_y] += da / na;
                        Sa[thread_y]   += da * (a - mu_a[thread_y]); 

                        PROPR_UNROLL
                        for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                            auto b = frag_b[1][thread_x];
                            auto db = b  - mu_b[thread_x];
                            mu_b[thread_x] += db / nb;
                            Sb[thread_x]   += db * (b - mu_b[thread_x]); 
                            accum[thread_y][thread_x] += da * (b - mu_b[thread_x]);
                            nb+=1;
                        }
                        na +=1;
                    }
                } while(tile_idx < K);
                
                const int c_block_row = a_tile_index;
                const int c_block_col = b_tile_index;
                const float denom = float(na + norm_type) -  1.0f;

                for (int i = 0; i < 4; ++i) {
                    float4 tmp0 = make_float4(
                        accum[i][0] / denom, accum[i][1]  / denom, 
                        accum[i][2] / denom, accum[i][3]  / denom
                    );

                    float4 tmp1 = make_float4(
                        accum[i][4] / denom, accum[i][5] / denom,
                        accum[i][6] / denom, accum[i][7] / denom
                    );

                    float4 tmp2 = make_float4(
                        accum[i+4][0] / denom, accum[i+4][1] / denom,
                        accum[i+4][2] / denom, accum[i+4][3] / denom
                    );

                    float4 tmp3 = make_float4(
                        accum[i+4][4] / denom, accum[i+4][5] / denom,
                        accum[i+4][6] / denom, accum[i+4][7] / denom
                    );

                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col,      M)]) = tmp0;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col + 64, M)]) = tmp1;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col,      M)]) = tmp2;
                    FETCH_FLOAT4(C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col + 64, M)]) = tmp3;
                }
            };

            __global__
            void 
            vlrRcpp(){

            };

            __global__
            void clrRcpp(
                    float * out,
                    offset_t out_stride,
                    float * __restrict__ x, 
                    offset_t x_stride,
                    size_t rows, size_t cols) {
                // TODO: investigate pipline sol (with split warps prod-cons)
                // TODO: f4
                const int col = blockDim.x * blockIdx.x + threadIdx.x;
                if ((size_t)col >= cols) return;

                float mean = 0.0;
                for (size_t r = 0; r < rows; ++r) {
                    float v = log2(x[r + col * x_stride]);
                    mean += (v - mean) / (r + 1);
                }
                for (size_t r = 0; r < rows; ++r) {
                    float v = x[r + col * x_stride];
                    out[r + col * out_stride] = (v - mean);
                }
            };

            __global__
            void alrRcpp(
                    const int ivar, 
                    float * out,
                    offset_t out_stride,
                    float * __restrict__ x, 
                    offset_t x_stride,
                    size_t rows, size_t cols) {
                // TODO: investigate pipline sol (with split warps prod-cons)
                // TODO: f4
                const int col = blockDim.x * blockIdx.x + threadIdx.x;
                if ((size_t)col >= cols) return;

                for (size_t r = 0; r < rows; ++r) {
                    out[r + col * out_stride] = log2(x[r + col * x_stride]) - log2(x[r + (ivar - 1) * x_stride]);
                }
                
            };

            __global__
            void symRcpp(){

            };

            __global__
            void phiRcpp(){

            };

            __global__
            void rhoRcpp(){

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

            __global__
            void linRcpp(){

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

                PROPR_UNROLL
                for (offset_t k = blockDim.x * blockIdx.x + threadIdx.x; k < total_pairs; k += gridDim.x * blockDim.x) {
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