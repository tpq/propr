#pragma once

#include <cub/cub.cuh>
#include <cuda_runtime.h>

#include <propr/data/traits.cuh>
#include <propr/utils/common/preprocessor.cuh>
#include <propr/internal/device/cuda/thread/mem_ops.cuh>



using namespace propr::cuda::internal;


template <typename Real, class Config>
__global__
void omega_kernel( 
    Real * __restrict__ A,
    Real * __restrict__ C,
    const int M,
    const int K) {

    using Wide = propr::cuda_wide_vector_t<Real>;
    constexpr int Lanes = propr::cuda_wide_lanes_v<Real>;

    static_assert((Config::BLK_K % Lanes)        == 0, "Config::BLK_K must be multiple of the selected vector width.");
    static_assert((Config::BLK_M % Config::TH_Y) == 0, "Config::BLK_M % Config::TH_Y == 0");
    static_assert((Config::BLK_M % Config::TH_X) == 0, "Config::BLK_M % Config::TH_X == 0");

    Real*B  = A;

    int bx = blockIdx.x;
    int by = blockIdx.y;
    // if (bx > by) return;

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    
    const int THREAD_X_PER_BLOCK = Config::BLK_M / Config::TH_X;
    const int THREAD_Y_PER_BLOCK = Config::BLK_M / Config::TH_Y;
    const int THREAD_NUM_PER_BLOCK = THREAD_X_PER_BLOCK * THREAD_Y_PER_BLOCK;

    const int tid = ty * THREAD_X_PER_BLOCK + tx;

    __shared__ Real As[2][Config::BLK_K][Config::BLK_M];
    __shared__ Real Bs[2][Config::BLK_K][Config::BLK_M];

    Real accum [Config::NUM_ACC][Config::TH_Y][Config::TH_X] = {Real(0)};

    Real frag_a[2][Config::TH_Y];
    Real frag_b[2][Config::TH_X];

    const int ldg_num_a = Config::BLK_M * Config::BLK_K / (THREAD_NUM_PER_BLOCK * Lanes);
    const int ldg_num_b = Config::BLK_K * Config::BLK_M / (THREAD_NUM_PER_BLOCK * Lanes);
    Real ldg_a_reg[Lanes*ldg_num_a];
    Real ldg_b_reg[Lanes*ldg_num_b];

    const int A_TILE_THREAD_PER_ROW = Config::BLK_K / Lanes;
    const int B_TILE_THREAD_PER_ROW = Config::BLK_M / Lanes;

    const int A_TILE_ROW_START = tid / A_TILE_THREAD_PER_ROW;
    const int B_TILE_ROW_START = tid / B_TILE_THREAD_PER_ROW;

    const int A_TILE_COL = tid % A_TILE_THREAD_PER_ROW * Lanes;
    const int B_TILE_COL = tid % B_TILE_THREAD_PER_ROW * Lanes;

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
        int ldg_index = i / A_TILE_ROW_STRIDE * Lanes;
        thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[ldg_index],
            thread::load<Config::LoadModifer, Wide>(&A[OFFSET( A_TILE_ROW_START + i,  A_TILE_COL, K )])
        );

        PROPR_UNROLL
        for (int lane = 0; lane < Lanes; ++lane) {
            As[0][A_TILE_COL + lane][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index + lane];
        }
    }

    PROPR_UNROLL
    for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
        int row_k = B_TILE_ROW_START + i;
        int col_n = B_TILE_COL;

        PROPR_UNROLL
        for (int lane = 0; lane < Lanes; ++lane) {
            Bs[0][row_k][col_n + lane] = B[OFFSET(col_n + lane, row_k, K)];
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

        if(tile_idx < K){
            PROPR_UNROLL
            for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                int ldg_index = i / A_TILE_ROW_STRIDE * Lanes;
                thread::store<Config::StoreModifer, Wide>(&ldg_a_reg[ldg_index],
                    thread::load<Config::LoadModifer, Wide>(&A[OFFSET(A_TILE_ROW_START + i,  A_TILE_COL + tile_idx, K )] )
                );
            }
            PROPR_UNROLL
            for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                int ldg_index = i / B_TILE_ROW_STRIDE * Lanes;
                int row_k = tile_idx + B_TILE_ROW_START + i;
                int col_n = B_TILE_COL;

                PROPR_UNROLL
                for (int lane = 0; lane < Lanes; ++lane) {
                    ldg_b_reg[ldg_index + lane] = B[OFFSET(col_n + lane, row_k, K)];
                }
            }
        }

        int load_stage_idx = write_stage_idx ^ 1;

        PROPR_UNROLL
        for(int j=0; j < Config::BLK_K - 1; ++j){
            PROPR_UNROLL
            for (int base = 0; base < Config::TH_Y; base += Lanes) {
                const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
                thread::store<Config::StoreModifer, Wide>(&frag_a[(j+1)%2][base],
                    thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx][(j+1)][a_tile_index + offset]));
            }

            PROPR_UNROLL
            for (int base = 0; base < Config::TH_X; base += Lanes) {
                const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
                thread::store<Config::StoreModifer, Wide>(&frag_b[(j+1)%2][base],
                    thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx][(j+1)][b_tile_index + offset]));
            }
            PROPR_UNROLL
            for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
                PROPR_UNROLL
                for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                    auto a = frag_a[j%2][thread_y];
                    auto b = frag_b[j%2][thread_x];
                    Config::template update<Real>(thread_x, thread_y, accum, a, b);
                }
            }
        }

        if(tile_idx < K){
            PROPR_UNROLL
            for ( int i = 0 ; i < Config::BLK_M ; i += A_TILE_ROW_STRIDE) {
                int ldg_index = i / A_TILE_ROW_STRIDE * Lanes;
                PROPR_UNROLL
                for (int lane = 0; lane < Lanes; ++lane) {
                    As[write_stage_idx][A_TILE_COL + lane][A_TILE_ROW_START + i] = ldg_a_reg[ldg_index + lane];
                }
            }
            PROPR_UNROLL
            for ( int i = 0 ; i < Config::BLK_K; i += B_TILE_ROW_STRIDE) {
                int ldg_index = i / B_TILE_ROW_STRIDE * Lanes;
                thread::store<Config::StoreModifer, Wide>(&Bs[write_stage_idx][B_TILE_ROW_START + i][B_TILE_COL],
                    thread::load<Config::LoadModifer, Wide>(&ldg_b_reg[ldg_index])
                );
            }
            __syncthreads();
            write_stage_idx ^= 1;
        }

        PROPR_UNROLL
        for (int base = 0; base < Config::TH_Y; base += Lanes) {
            const int offset = (base < Config::TH_Y / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_Y / 2);
            thread::store<Config::StoreModifer, Wide>(&frag_a[0][base],
                thread::load<Config::LoadModifer, Wide>(&As[load_stage_idx^1][0][a_tile_index + offset]));
        }

        PROPR_UNROLL
        for (int base = 0; base < Config::TH_X; base += Lanes) {
            const int offset = (base < Config::TH_X / 2) ? base : Config::BLK_M / 2 + (base - Config::TH_X / 2);
            thread::store<Config::StoreModifer, Wide>(&frag_b[0][base],
                thread::load<Config::LoadModifer, Wide>(&Bs[load_stage_idx^1][0][b_tile_index + offset]));
        }

        PROPR_UNROLL
        for (int thread_y = 0; thread_y < Config::TH_Y; ++thread_y) {
            PROPR_UNROLL
            for (int thread_x = 0; thread_x < Config::TH_X; ++thread_x) {
                auto a = frag_a[1][thread_y];
                auto b = frag_b[1][thread_x];
                Config::template update<Real>(thread_x, thread_y, accum, a, b);
            }
        }
    } while(tile_idx< K);
    
    const int c_block_row = a_tile_index;
    const int c_block_col = b_tile_index;

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
            thread::store<Config::StoreModifer, Real>(
                &C[OFFSET(row, col, M)],
                Config::template finalize<Real>(thread_x, thread_y, accum)
            );
        }
    }

}
