#pragma once

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <propr/utils/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>



using namespace propr::cuda::internal;


template <class Config> 
__global__
void omega_kernel( 
    float * __restrict__ A,
    float * __restrict__ C, 
    const int M,
    const int K) {

    static_assert((Config::BLK_K % 4)             == 0, "Config::BLK_K must be multiple of 4 (float4 gmem loads).");
    static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
    static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

    float*B  = A;

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

    float accum [Config::NUM_ACC][Config::TH_Y][Config::TH_X] = {0.0f};

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
    const int a_tile_index =  (warp_id / 2) * 16 + (lane_id / 8) * 4;
    const int b_tile_index =  (warp_id % 2) * 32 + (lane_id % 8) * 4;
    
    PROPR_UNROLL
    for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
        int ldg_index = i / A_TILE_ROW_STRIDE * 4;
        thread::store<Config::StoreModifer,float4>(&ldg_a_reg[ldg_index],  
            thread::load<Config::LoadModifer,float4>(&A[OFFSET( A_TILE_ROW_START + i,  A_TILE_COL, K )])
        );

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
    
    thread::store<Config::StoreModifer,float4>(&frag_a[0][0], 
        thread::load<Config::LoadModifer,float4>(&As[0][0][a_tile_index]));
    thread::store<Config::StoreModifer,float4>(&frag_a[0][4], 
        thread::load<Config::LoadModifer,float4>(&As[0][0][a_tile_index + 64]));
    
    thread::store<Config::StoreModifer,float4>(&frag_b[0][0], 
        thread::load<Config::LoadModifer,float4>(&Bs[0][0][b_tile_index]));
    thread::store<Config::StoreModifer,float4>(&frag_b[0][4],
         thread::load<Config::LoadModifer,float4>(&Bs[0][0][b_tile_index + 64]));
    
    int write_stage_idx = 1;
    int tile_idx = 0;
    do {
        tile_idx += Config::BLK_K;

        if(tile_idx < K){
            PROPR_UNROLL
            for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                thread::store<Config::StoreModifer,float4>(&ldg_a_reg[ldg_index], 
                    thread::load<Config::LoadModifer,float4>(&A[OFFSET(A_TILE_ROW_START + i,  A_TILE_COL + tile_idx, K )] )
                );
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
            thread::store<Config::StoreModifer,float4>(&frag_a[(j+1)%2][0],  
                thread::load<Config::LoadModifer,float4>(&As[load_stage_idx][(j+1)][a_tile_index]));
            thread::store<Config::StoreModifer,float4>(&frag_a[(j+1)%2][4],  
                thread::load<Config::LoadModifer,float4>(&As[load_stage_idx][(j+1)][a_tile_index + 64])
            );

            thread::store<Config::StoreModifer,float4>(&frag_b[(j+1)%2][0],  
                thread::load<Config::LoadModifer,float4>(&Bs[load_stage_idx][(j+1)][b_tile_index])
            );
            thread::store<Config::StoreModifer,float4>(&frag_b[(j+1)%2][4],  
                thread::load<Config::LoadModifer,float4>(&Bs[load_stage_idx][(j+1)][b_tile_index + 64])
            );
            PROPR_UNROLL
            for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                PROPR_UNROLL
                for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                    auto a = frag_a[j%2][thread_y];
                    auto b = frag_b[j%2][thread_x];
                    Config::update(thread_x, thread_y, accum, a, b);
                }
            }
        }

        if(tile_idx < K){
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
                thread::store<Config::StoreModifer,float4>(&Bs[write_stage_idx][B_TILE_ROW_START + i][B_TILE_COL], 
                    thread::load<Config::LoadModifer,float4>(&ldg_b_reg[ldg_index])
                );
            }
            __syncthreads();
            write_stage_idx ^= 1;
        }

        thread::store<Config::StoreModifer,float4>(&frag_a[0][0],
            thread::load<Config::LoadModifer,float4>(&As[load_stage_idx^1][0][a_tile_index]));
        thread::store<Config::StoreModifer,float4>(&frag_a[0][4],
            thread::load<Config::LoadModifer,float4>(&As[load_stage_idx^1][0][a_tile_index + 64]));
        
        thread::store<Config::StoreModifer,float4>(&frag_b[0][0],
            thread::load<Config::LoadModifer,float4>(&Bs[load_stage_idx^1][0][b_tile_index]));
        thread::store<Config::StoreModifer,float4>(&frag_b[0][4],
            thread::load<Config::LoadModifer,float4>(&Bs[load_stage_idx^1][0][b_tile_index + 64]));

        PROPR_UNROLL
        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
            PROPR_UNROLL
            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                auto a = frag_a[1][thread_y];
                auto b = frag_b[1][thread_x];
                Config::update(thread_x, thread_y, accum, a, b);
            }
        }
    } while(tile_idx< K);
    
    const int c_block_row = a_tile_index;
    const int c_block_col = b_tile_index;

    for (int i = 0; i < 4; ++i) {
        float4 tmp0 = make_float4(
            Config::finalize(i,0, accum),
            Config::finalize(i,1, accum),
            Config::finalize(i,2, accum),
            Config::finalize(i,3, accum)
        );

        float4 tmp1 = make_float4(
            Config::finalize(i, 4, accum),
            Config::finalize(i, 5, accum),
            Config::finalize(i, 6, accum),
            Config::finalize(i, 7, accum)
        );

        float4 tmp2 = make_float4(
            Config::finalize(i + 4, 0, accum),
            Config::finalize(i + 4, 1, accum),
            Config::finalize(i + 4, 2, accum),
            Config::finalize(i + 4, 3, accum)
        );

        float4 tmp3 = make_float4(
            Config::finalize(i + 4, 4, accum),
            Config::finalize(i + 4, 5, accum),
            Config::finalize(i + 4, 6, accum),
            Config::finalize(i + 4, 7, accum)
        );

        thread::store<Config::StoreModifer,float4>(&C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col,      M)], thread::load<Config::LoadModifer,float4>(&tmp0));
        thread::store<Config::StoreModifer,float4>(&C[OFFSET(Config::BLK_M * by + c_block_row + i,      Config::BLK_M * bx + c_block_col + 64, M)], thread::load<Config::LoadModifer,float4>(&tmp1));
        thread::store<Config::StoreModifer,float4>(&C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col,      M)], thread::load<Config::LoadModifer,float4>(&tmp2));
        thread::store<Config::StoreModifer,float4>(&C[OFFSET(Config::BLK_M * by + c_block_row + 64 + i, Config::BLK_M * bx + c_block_col + 64, M)], thread::load<Config::LoadModifer,float4>(&tmp3));
    }

}