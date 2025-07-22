#pragma once

#include <cute/tensor.hpp>
#include <propr/kernels/cuda/detail/cutlass/omega/thread/statistics.cuh>

namespace propr {
    namespace kernels {
        namespace cutlass {
            template <typename Config, typename T=cute::half_t>
            __global__ 
            void omega_kernel(const T *Aptr, T *Dptr, int m, int k) {
                using namespace cute;

                int idx = threadIdx.x;
                int ix  = blockIdx.x;
                int iy  = blockIdx.y;
                if (ix > iy) return;

                Tensor A  = make_tensor(make_gmem_ptr(Aptr), make_shape(m, k), make_stride(k, Int<1>{}));
                Tensor B  = make_tensor(make_gmem_ptr(Aptr), make_shape(m, k), make_stride(k, Int<1>{})); // non-owing copy of A
                Tensor D  = make_tensor(make_gmem_ptr(Dptr), make_shape(m, m), make_stride(m, Int<1>{}));
                Tensor DT = make_tensor(make_gmem_ptr(Dptr), make_shape(m, m), make_stride(Int<1>{}, m));  // non-owing copy of D transpose

                Tensor gA  = local_tile(A , make_tile(Int<Config::BLK_M>{}, Int<Config::BLK_K>{}), make_coord(iy,  _)); // (BLK_M, BLK_K, num_tile_k)
                Tensor gB  = local_tile(B , make_tile(Int<Config::BLK_M>{}, Int<Config::BLK_K>{}), make_coord(ix,  _)); // (BLK_M, BLK_K, num_tile_k)
                Tensor gD  = local_tile(D , make_tile(Int<Config::BLK_M>{}, Int<Config::BLK_M>{}), make_coord(iy, ix)); // (BLK_M, BLK_M)
                Tensor gDT = local_tile(DT, make_tile(Int<Config::BLK_M>{}, Int<Config::BLK_M>{}), make_coord(iy, ix)); // (BLK_M, BLK_M)

                // Compute tile residues for predication
                auto m_max_coord = m - size<0>(gA) * iy; 
                auto n_max_coord = m - size<0>(gB) * ix; 
                auto k_residue   = k - size<1>(gA) * (size<2>(gA) - 1);
                auto residue_mnk = make_tuple(m_max_coord, n_max_coord, k_residue);

                // Shared memory buffers
                // use char type to avoid issues with different template param T
                extern __shared__ char shared_memory[];
                Config::SharedStorage &smem = *reinterpret_cast<Config::SharedStorage *>(shared_memory);
                auto sA = make_tensor(make_smem_ptr(smem.A.begin()), Config::SmemLayoutA{}); // (BLK_M,BLK_K,PIPE)
                auto sB = make_tensor(make_smem_ptr(smem.B.begin()), Config::SmemLayoutB{});

                // register, use tiled_mma to partition register A/B/C
                typename Config::TiledMMA tiled_mma;
                auto thr_mma = tiled_mma.get_slice(threadIdx.x);
                auto tCgD    = thr_mma.partition_C(gD); 
                auto tCgDT   = thr_mma.partition_C(gDT); // partition for D transposed

                auto tCrA   = thr_mma.partition_fragment_A(gA(_, _, 0)); // (MMA, MMA_M, MMA_K)
                auto tCrB   = thr_mma.partition_fragment_B(gB(_, _, 0)); // (MMA, MMA_N, MMA_K)

                thread::OmegaHarmonic<Config> thr_omg(gD);

                // auto tCrD_n = thr_mma.partition_fragment_C(gD);          // (MMA, MMA_M, MMA_N)
                // auto tCrD_s = thr_mma.partition_fragment_C(gD);
                // clear(tCrD_n); clear(tCrD_s);

                // from global memory to shared memory
                typename Config::GmemCopyA g2s_tiled_copy_a;
                typename Config::GmemCopyB g2s_tiled_copy_b;
                auto g2s_thr_copy_a = g2s_tiled_copy_a.get_slice(idx);
                auto g2s_thr_copy_b = g2s_tiled_copy_b.get_slice(idx);
                auto tAgA_copy = g2s_thr_copy_a.partition_S(gA); // (CPY, CPY_M, CPY_K, k)
                auto tAsA_copy = g2s_thr_copy_a.partition_D(sA); // (CPY, CPY_M, CPY_K)
                auto tBgB_copy = g2s_thr_copy_b.partition_S(gB); // (CPY, CPY_N, CPY_K, k)
                auto tBsB_copy = g2s_thr_copy_b.partition_D(sB); // (CPY, CPY_N, CPY_K)

                // from shared memory to register, use tiled_mma to generate tiled_copy
                typename Config::SmemCopyA s2r_tiled_copy_a;

                auto s2r_thr_copy_a   = s2r_tiled_copy_a.get_slice(idx);
                auto tAsA = s2r_thr_copy_a.partition_S(sA);     // (CPY, CPY_M, CPY_K)
                auto tCrA_view = s2r_thr_copy_a.retile_D(tCrA); // (CPY, CPY_M, CPY_K)

                typename Config::SmemCopyB s2r_tiled_copy_b;
                auto s2r_thr_copy_b   = s2r_tiled_copy_b.get_slice(idx);
                auto tBsB = s2r_thr_copy_b.partition_S(sB);     // (CPY, CPY_N, CPY_K)
                auto tCrB_view = s2r_thr_copy_b.retile_D(tCrB); // (CPY, CPY_N, CPY_K)

                //
                // PREDICATES
                //
                Tensor tApA = make_tensor<bool>(make_shape(size<1>(tAsA_copy), size<2>(tAsA_copy)), Stride<_1, _0>{});
                Tensor tBpB = make_tensor<bool>(make_shape(size<1>(tBsB_copy), size<2>(tBsB_copy)), Stride<_1, _0>{});

                Tensor cA = make_identity_tensor(make_shape(size<0>(sA), size<1>(sA)));
                Tensor cB = make_identity_tensor(make_shape(size<0>(sB), size<1>(sB)));

                Tensor tAcA = g2s_thr_copy_a.partition_S(cA);
                Tensor tBcB = g2s_thr_copy_b.partition_S(cB);

                CUTE_UNROLL
                for (int i = 0; i < size<0>(tApA); ++i) {
                    tApA(i, 0) = get<0>(tAcA(0, i, 0)) < m_max_coord; // blk_m coord < m_max_coord
                }
            
                CUTE_UNROLL
                for (int i = 0; i < size<0>(tBpB); ++i) {
                    tBpB(i, 0) = get<0>(tBcB(0, i, 0)) < n_max_coord; // blk_n coord < n_max_coord
                }

                // Clear the smem tiles to account for predicated off loads
                clear(tAsA_copy);
                clear(tBsB_copy);

                Tensor tAgAk = tAgA_copy(_, _, _, size<2>(gA) - 1);
                CUTE_UNROLL
                for (int i = 0; i < size<2>(tAsA_copy); ++i) {
                    if (get<1>(tAcA(0, 0, i)) < get<2>(residue_mnk)) { 
                        // blk_k coord < residue_k (gA shifted)
                        cute::copy_if(g2s_tiled_copy_a, tApA(_, i), tAgAk(_, _, i), tAsA_copy(_, _, i));
                    }
                }

                Tensor tBgBLK_K = tBgB_copy(_, _, _, size<2>(gB) - 1);
                CUTE_UNROLL
                for (int i = 0; i < size<2>(tBsB_copy); ++i) {
                    if (get<1>(tBcB(0, 0, i)) < get<2>(residue_mnk)) { 
                        // blk_k coord < residue_k (gA shifted)
                        cute::copy_if(g2s_tiled_copy_b, tBpB(_, i), tBgBLK_K(_, _, i), tBsB_copy(_, _, i));
                    }
                }
                cp_async_fence();

                int ntile = (k + Config::BLK_K - 1) / Config::BLK_K;
                CUTE_UNROLL
                for (int itile = 0; itile < ntile; ++itile) {
                    if (itile >= 1) {
                        cute::copy_if(g2s_tiled_copy_a, tApA, tAgA_copy(_, _, _, itile - 1), tAsA_copy(_, _, _));
                        cute::copy_if(g2s_tiled_copy_b, tBpB, tBgB_copy(_, _, _, itile - 1), tBsB_copy(_, _, _));
                        cp_async_fence();
                    }

                    cp_async_wait<0>();
                    __syncthreads();

                    int nk = size<2>(tCrA);
                    CUTE_UNROLL
                    for (int ik = 0; ik < nk; ++ik) {
                        cute::copy(s2r_tiled_copy_a, tAsA(_, _, ik), tCrA_view(_, _, ik)); // copy  (CPY, CPY_M), sync
                        cute::copy(s2r_tiled_copy_b, tBsB(_, _, ik), tCrB_view(_, _, ik)); // copy  (CPY, CPY_N)
                        auto a_view = tCrA(_, _, ik); auto b_view = tCrB(_, _, ik);
                        thr_omg(a_view, b_view);
                    }
                } 

                __syncthreads();

                Tensor cD = make_identity_tensor(make_shape(size<0>(gD), size<1>(gD)));
                Tensor tCcD = thr_mma.partition_C(cD); 
                typename Config::NumericConverter converter;
                CUTE_UNROLL
                for (int i = 0; i < thr_omg.size(); ++i) {
                    if (elem_less(tCcD(i), make_coord(m_max_coord, n_max_coord))) {
                        auto result = thr_omg.get(i);
                        tCgD(i) = converter(result);
                        if (ix < iy) tCgDT(i) = converter(result);
                    }
                }

            }

        }
    }
}