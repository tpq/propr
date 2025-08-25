#pragma once

#include <cuda_runtime.h>
#include <propr/utils/preprocessor.cuh>

#define FETCH_FLOAT4(pointer) (reinterpret_cast<float4*>(&(pointer))[0])
#define OFFSET(row, col, ld) ((row) * (ld) + (col))


template <
    const int BLK_M,
    const int BLK_K,
    const int TH_Y,
    const int TH_X> 
__global__
void omega_kernel( 
    float * __restrict__ A,
    float * __restrict__ C, 
    const int M,
    const int K) {

    float*B  = A;

    int bx = blockIdx.x;
    int by = blockIdx.y;
    // if (bx > by) return;

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    
    const int THREAD_X_PER_BLOCK = BLK_M / TH_X;
    const int THREAD_Y_PER_BLOCK = BLK_M / TH_Y;
    const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;

    const int tid = ty * THREAD_X_PER_BLOCK + tx;

    __shared__ float As[2][BLK_K][BLK_M];
    __shared__ float Bs[2][BLK_K][BLK_M];

    float accum [TH_Y][TH_X] = {0.0f};
    float accum2[TH_Y][TH_X] = {0.0f};

    float frag_a[2][TH_Y];
    float frag_b[2][TH_X];

    const int ldg_num_a = BLK_M * BLK_K / (THREAD_NUM_PER_BLOCK * 4);
    const int ldg_num_b = BLK_K * BLK_M / (THREAD_NUM_PER_BLOCK * 4);
    float ldg_a_reg[4*ldg_num_a];
    float ldg_b_reg[4*ldg_num_b];

    const int A_TILE_THREAD_PER_ROW = BLK_K / 4;
    const int B_TILE_THREAD_PER_ROW = BLK_M / 4;

    const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
    const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

    const int A_TILE_COL = tid % A_TILE_THREAD_PER_ROW * 4;
    const int B_TILE_COL = tid % B_TILE_THREAD_PER_ROW * 4;

    const int A_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / A_TILE_THREAD_PER_ROW;
    const int B_TILE_ROW_STRIDE = THREAD_NUM_PER_BLOCK / B_TILE_THREAD_PER_ROW;

    A = &A[(BLK_M * by) * K];
    B = &B[(BLK_M * bx) * K];

    const int warp_id = tid / 32;
    const int lane_id = tid % 32;
    const int a_tile_index =  warp_id/2*16 + lane_id/8*4;
    const int b_tile_index =  warp_id%2*32 + lane_id%8*4;
    
    PROPR_UNROLL
    for ( int i = 0 ; i < BLK_M ; i += A_TILE_ROW_STRIDE) {
        int ldg_index = i / A_TILE_ROW_STRIDE * 4;
        FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(
            A_TILE_ROW_START + i, // row
            A_TILE_COL, // col
            K )]);
        As[0][A_TILE_COL  ][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index  ];
        As[0][A_TILE_COL+1][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+1];
        As[0][A_TILE_COL+2][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+2];
        As[0][A_TILE_COL+3][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+3];
    }

    PROPR_UNROLL
    for ( int i = 0 ; i < BLK_K; i += B_TILE_ROW_STRIDE) {
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
        tile_idx += BLK_K;

        if(tile_idx < K){
            PROPR_UNROLL
            for ( int i = 0 ; i < BLK_M ; i += A_TILE_ROW_STRIDE) {
                int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                FETCH_FLOAT4(ldg_a_reg[ldg_index]) = FETCH_FLOAT4(A[OFFSET(
                    A_TILE_ROW_START + i, // row
                    A_TILE_COL + tile_idx, // col
                    K )]);
            }
            PROPR_UNROLL
            for ( int i = 0 ; i < BLK_K; i += B_TILE_ROW_STRIDE) {
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
        for(int j=0; j<BLK_K - 1; ++j){
            FETCH_FLOAT4(frag_a[(j+1)%2][0]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index]);
            FETCH_FLOAT4(frag_a[(j+1)%2][4]) = FETCH_FLOAT4(As[load_stage_idx][(j+1)][a_tile_index + 64]);
            FETCH_FLOAT4(frag_b[(j+1)%2][0]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index]);
            FETCH_FLOAT4(frag_b[(j+1)%2][4]) = FETCH_FLOAT4(Bs[load_stage_idx][(j+1)][b_tile_index + 64]);
            PROPR_UNROLL
            for (int thread_y = 0; thread_y < TH_Y; ++thread_y) {
                PROPR_UNROLL
                for (int thread_x = 0; thread_x < TH_X; ++thread_x) {
                    auto a = frag_a[j%2][thread_y];
                    auto b = frag_b[j%2][thread_x];
                    auto n = 2.0f * (a / ( a + b + FLT_EPSILON)) * b;
                    accum[thread_y][thread_x]  += n;
                    accum2[thread_y][thread_x] += n * n;
                }
            }
        }

        if(tile_idx < K){
            PROPR_UNROLL
            for ( int i = 0 ; i < BLK_M ; i += A_TILE_ROW_STRIDE) {
                int ldg_index = i / A_TILE_ROW_STRIDE * 4;
                As[write_stage_idx][A_TILE_COL  ][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index];
                As[write_stage_idx][A_TILE_COL+1][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+1];
                As[write_stage_idx][A_TILE_COL+2][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+2];
                As[write_stage_idx][A_TILE_COL+3][A_TILE_ROW_START + i]=ldg_a_reg[ldg_index+3];
            }
            PROPR_UNROLL
            for ( int i = 0 ; i < BLK_K; i += B_TILE_ROW_STRIDE) {
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
        for (int thread_y = 0; thread_y < TH_Y; ++thread_y) {
            PROPR_UNROLL
            for (int thread_x = 0; thread_x < TH_X; ++thread_x) {
                auto a = frag_a[1][thread_y];
                auto b = frag_b[1][thread_x];
                auto n = 2.0f * (a / ( a + b + FLT_EPSILON)) * b;
                accum [thread_y][thread_x] += n;
                accum2[thread_y][thread_x] += n * n;
            }
        }
    } while(tile_idx< K);
    
    const int c_block_row = a_tile_index;
    const int c_block_col = b_tile_index;

    for (int i = 0; i < 4; ++i) {
        float4 tmp0 = make_float4(
            accum[i][0] - accum2[i][0] / accum[i][0],
            accum[i][1] - accum2[i][1] / accum[i][1],
            accum[i][2] - accum2[i][2] / accum[i][2],
            accum[i][3] - accum2[i][3] / accum[i][3]
        );

        float4 tmp1 = make_float4(
            accum[i][4] - accum2[i][4] / accum[i][4],
            accum[i][5] - accum2[i][5] / accum[i][5],
            accum[i][6] - accum2[i][6] / accum[i][6],
            accum[i][7] - accum2[i][7] / accum[i][7]
        );

        float4 tmp2 = make_float4(
            accum[i+4][0] - accum2[i+4][0] / accum[i+4][0],
            accum[i+4][1] - accum2[i+4][1] / accum[i+4][1],
            accum[i+4][2] - accum2[i+4][2] / accum[i+4][2],
            accum[i+4][3] - accum2[i+4][3] / accum[i+4][3]
        );

        float4 tmp3 = make_float4(
            accum[i+4][4] - accum2[i+4][4] / accum[i+4][4],
            accum[i+4][5] - accum2[i+4][5] / accum[i+4][5],
            accum[i+4][6] - accum2[i+4][6] / accum[i+4][6],
            accum[i+4][7] - accum2[i+4][7] / accum[i+4][7]
        );

        FETCH_FLOAT4(C[OFFSET(BLK_M * by + c_block_row + i,           BLK_M * bx + c_block_col,          M)]) = tmp0;
        FETCH_FLOAT4(C[OFFSET(BLK_M * by + c_block_row + i,           BLK_M * bx + c_block_col + 64,     M)]) = tmp1;
        FETCH_FLOAT4(C[OFFSET(BLK_M * by + c_block_row + 64 + i,      BLK_M * bx + c_block_col,          M)]) = tmp2;
        FETCH_FLOAT4(C[OFFSET(BLK_M * by + c_block_row + 64 + i,      BLK_M * bx + c_block_col + 64,     M)]) = tmp3;
    }

}