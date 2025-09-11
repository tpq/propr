#pragma once

#include <cute/tensor.hpp>
#include <propr/utils/constants.h>

namespace propr {
    namespace kernels {
        namespace cutlass_impl {
            template <typename Config, typename Tin, typename Tout>
            __global__  
            // __launch_bounds__(int(decltype(size(Config::TiledMMA{}))::value))
            void omega_kernel(int m, int k, const Tin *Aptr, Tout *Dptr) {
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
                typename Config::SharedStorage &smem = *reinterpret_cast<typename Config::SharedStorage *>(shared_memory);
                typename Config::SmemLayoutA sA_layout;
                typename Config::SmemLayoutB sB_layout;

                auto sA = make_tensor(make_smem_ptr(smem.A.begin()), sA_layout); // (BLK_M,BLK_K,PIPE)
                auto sB = make_tensor(make_smem_ptr(smem.B.begin()), sB_layout);

                // register, use tiled_mma to partition register A/B/C
                typename Config::TiledMMA  tiled_mma;
                typename Config::TiledMMA2 tiled_mma2;

                auto thr_mma = tiled_mma.get_slice(threadIdx.x);
                auto tCgD    = thr_mma.partition_C(gD); 
                auto tCgDT   = thr_mma.partition_C(gDT); // partition for D transposed

                auto tCrA   = thr_mma.partition_fragment_A(gA(_, _, 0)); // (MMA, MMA_M, MMA_K)
                auto tCrB   = thr_mma.partition_fragment_B(gB(_, _, 0)); // (MMA, MMA_N, MMA_K)

                auto tCrD_n  = thr_mma.partition_fragment_C(gD);
                auto tCrD_s  = thr_mma.partition_fragment_C(gD);
                auto dummy   = thr_mma.partition_fragment_C(gD);

                clear(tCrD_n); clear(tCrD_s); clear(dummy);

                // from global memory to shared memory
                typename Config::GmemCopyA g2s_tiled_copy_a;
                auto g2s_thr_copy_a = g2s_tiled_copy_a.get_slice(idx);
                auto tAgA_copy      = g2s_thr_copy_a.partition_S(gA); // (CPY, CPY_M, CPY_K, k)
                auto tAsA_copy      = g2s_thr_copy_a.partition_D(sA); // (CPY, CPY_M, CPY_K)

                typename Config::GmemCopyB g2s_tiled_copy_b;
                auto g2s_thr_copy_b = g2s_tiled_copy_b.get_slice(idx);
                auto tBgB_copy      = g2s_thr_copy_b.partition_S(gB); // (CPY, CPY_N, CPY_K, k)
                auto tBsB_copy      = g2s_thr_copy_b.partition_D(sB); // (CPY, CPY_N, CPY_K)

                // from shared memory to register, use tiled_mma to generate tiled_copy
                typename Config::SmemCopyA s2r_tiled_copy_a;
                auto s2r_thr_copy_a   = s2r_tiled_copy_a.get_slice(idx);
                auto tAsA             = s2r_thr_copy_a.partition_S(sA);     // (CPY, CPY_M, CPY_K)
                auto tCrA_view        = s2r_thr_copy_a.retile_D(tCrA);      // (CPY, CPY_M, CPY_K)

                typename Config::SmemCopyB s2r_tiled_copy_b;
                auto s2r_thr_copy_b = s2r_tiled_copy_b.get_slice(idx);
                auto tBsB           = s2r_thr_copy_b.partition_S(sB);     // (CPY, CPY_N, CPY_K)
                auto tCrB_view      = s2r_thr_copy_b.retile_D(tCrB);      // (CPY, CPY_N, CPY_K)

                //
                // PREDICATES
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
                clear(tAsA_copy); clear(tBsB_copy);

                Tensor tAgAk = tAgA_copy(_, _, _, size<2>(gA) - 1);
                CUTE_UNROLL
                for (int i = 0; i < size<2>(tAsA_copy); ++i) {
                    if (get<1>(tAcA(0, 0, i)) < get<2>(residue_mnk)) { // blk_k coord < residue_k (gA shifted)
                        cute::copy_if(g2s_tiled_copy_a, tApA(_, i), tAgAk(_, _, i), tAsA_copy(_, _, i));
                    }
                }

                Tensor tBgBk = tBgB_copy(_, _, _, size<2>(gB) - 1);
                CUTE_UNROLL
                for (int i = 0; i < size<2>(tBsB_copy); ++i) {
                    if (get<1>(tBcB(0, 0, i)) < get<2>(residue_mnk)) {  // blk_k coord < residue_k (gA shifted)
                        cute::copy_if(g2s_tiled_copy_b, tBpB(_, i), tBgBk(_, _, i), tBsB_copy(_, _, i));
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

                    // int nk = size<2>(tCrA);
                    // CUTE_UNROLL
                    // for (int ik = 0; ik < nk; ++ik) {
                    //     cute::copy(s2r_tiled_copy_a, tAsA(_, _, ik), tCrA_view(_, _, ik)); // copy  (CPY, CPY_M), sync
                    //     cute::copy(s2r_tiled_copy_b, tBsB(_, _, ik), tCrB_view(_, _, ik)); // copy  (CPY, CPY_N)

                    //     cute::gemm(tiled_mma , acc_1, tCrA(_, _, ik), tCrB(_, _, ik), acc_1);
                    //     cute::gemm(tiled_mma2, acc_2, tCrA(_, _, ik), tCrB(_, _, ik), acc_2);
                        
                    //     // printf("\n");
                    // }
                    // // printf("\n\n\n");

                    int nk = size<2>(tCrA);
                    CUTE_UNROLL
                    for (int ik = 0; ik < nk; ++ik) {
                        cute::copy(s2r_tiled_copy_a, tAsA(_, _, ik), tCrA_view(_, _, ik)); // copy  (CPY, CPY_M), sync
                        for(int ki=0; ki < 8 ; ki++) {
                            int warp_id  = threadIdx.x / PROPR_WARP_SIZE;
                            int lane_id  = threadIdx.x % PROPR_WARP_SIZE;
                            int remapped_indx =  warp_id * PROPR_WARP_SIZE + (ki* 4 + (lane_id % 4));
                            // printf(" threadidx.x = %d remapped = %d\n",threadIdx.x, remapped_indx);
                            auto s2r_thr_copy_b = s2r_tiled_copy_b.get_slice(idx);
                            // auto s2r_thr_copy_b = s2r_tiled_copy_b.get_slice(threadIdx.x);
                            auto tBsB           = s2r_thr_copy_b.partition_S(sB);     // (CPY, CPY_N, CPY_K)
                            auto tCrB_view      = s2r_thr_copy_b.retile_D(tCrB);
                             
                            cute::copy(s2r_tiled_copy_b, tBsB(_, _, ik), tCrB_view(_, _, ik)); // copy  (CPY, CPY_N)

                            auto& acc_1 =  (lane_id % 4) == (ki % 4) ? tCrD_n : dummy;
                            auto& acc_2 =  (lane_id % 4) == (ki % 4) ? tCrD_s : dummy;

                            if (threadIdx.x == 1) print_tensor(tCrB(_, _, ik));
                            __syncthreads();

                            cute::gemm(tiled_mma , acc_1, tCrA(_, _, ik), tCrB(_, _, ik), acc_1);
                            cute::gemm(tiled_mma2, acc_2, tCrA(_, _, ik), tCrB(_, _, ik), acc_2);

                            // if (threadIdx.x == 0|| threadIdx.x == 1 || threadIdx.x == 2 || threadIdx.x == 3) 
                            //     printf("%d[%d] %f %f \n",(lane_id % 4) == (ki % 4), threadIdx.x, tCrD_n(0), tCrD_s(0));
                        }
                        // printf("\n");
                    }
                    // printf("\n\n\n");
                    
                }

                __syncthreads();

                Tensor cD = make_identity_tensor(make_shape(size<0>(gD), size<1>(gD)));
                Tensor tCcD = thr_mma.partition_C(cD); 
                CUTE_UNROLL
                for (int i = 0; i < cute::size(tCrD_s); ++i) {
                    if (elem_less(tCcD(i), make_coord(m_max_coord, n_max_coord))) {
                        auto n_val = tCrD_n(i), s_val = tCrD_s(i);
                        if (threadIdx.x == 4) printf("[%d] %f %f \n",threadIdx.x, n_val, s_val);
                        auto result = n_val;
                        // auto result = s_val;
                        // auto result = (1.0f - (n_val == 0.0f)) * (n_val - s_val / (n_val + (n_val == 0.0f)));
                        // auto result = tCrD_n(i) - tCrD_s(i);
                        tCgD(i) = result;
                        if (ix < iy) tCgDT(i) = result;
                    }
                }

                printAgBg();
            }
        }
    }
}


// void 
// propr::dispatch::cuda::dof_global(NumericVector& out, const NumericMatrix& W, propr_context context) {
//     using Config = kernels::cutlass_impl::OmegaConfig;
//     NumericMatrix W2(W);
//     // printRMatrix(W2);
//     // for(int i=0; i < W2.nrow(); i++){
//     //   for(int j=0; j < W2.ncol(); j++){
//     //     W2(i,j) = float(335.712);
//     //   } 
//     // }

//     for(int i=0; i < W2.nrow(); i++){
//       for(int j=0; j < W2.ncol(); j++){
//         W2(i,j) = float(0);
//       } 
//     }

//     for(int i=0; i < W2.nrow(); i++){
//       for(int j=0; j < W2.ncol(); j++){
//         W2(i,j) = float(j) * W2.nrow() + i + 1;
//       }
//     }
//     printRMatrix(W2);

//     int t = 32;
//     auto Wl = pad_matrix(W2, 0, ((W2.nrow() + t - 1)/t)*t - W2.nrow(), 0, ((W2.ncol() + t - 1)/t)*t - W2.ncol());
//     int nfeats  = Wl.ncol();
//     int samples = Wl.nrow();
//     // printRMatrix(Wl);

//     //  for(int i=0; i < Wl.nrow(); i++){
//     //   for(int j=0; j < Wl.ncol(); j++){
//     //     Wl(i,j) = float(i) + float(j)*Wl.nrow();
//     //   }
//     // }
//     // printRMatrix(Wl);
    
//     std::cout << "[" << W.nrow() << "," << W.ncol() << "]" << std::endl;

//     const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
//     // CHECK_VECTOR_SIZE(out, llt);
//     int alignment = 1; //16
//     offset_t W_stride;
//     auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP>(Wl, W_stride, alignment);
//     float* d_out = nullptr;
//     offset_t dout_stride = ((nfeats  + alignment - 1) / alignment) * alignment;
//     CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));

//     const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
//     const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

//     dim3 block(cute::size(Config::TiledMMA{}));
//     dim3 grid(BX, BY);
//     static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
//     static constexpr int kShmSize = shm_size_AB * sizeof(__half);
//     const auto fptr = propr::kernels::cutlass_impl::omega_kernel<Config, cute::half_t, float>;

//     cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, kShmSize);
//     cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);

//     std::cout <<  "M: "<< nfeats << " K: " << samples << std::endl;
//     std::cout <<  "fptr<<<(" << grid.x << ","<< grid.y <<")"<< ",(" << block.x << "," << block.y <<")>>>" << std::endl;
//     fptr<<<grid, block, kShmSize, context.stream>>>(nfeats, samples, d_W, d_out);
    
//     CUDA_CHECK(cudaStreamSynchronize(context.stream));
//     auto h_full= new std::vector<float> (nfeats * dout_stride);
//     CUDA_CHECK(cudaMemcpy(
//         h_full->data(),
//         d_out,
//         nfeats * dout_stride * sizeof(float),
//         cudaMemcpyDeviceToHost
//     ));

//     size_t counter = 0;
//     double* out_ptr = REAL(out);
//     std::cout << "[GPU]: \n";
//     for (int i = 0; i < nfeats; ++i) {
//         for (int j = 0; j < nfeats; ++j) {
//             float v = h_full->at(size_t(i) * dout_stride + j);
//             // out_ptr[counter++] = static_cast<double>(v);
//             std::cout << v << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;

//     CUDA_CHECK(cudaFree(d_W));
//     CUDA_CHECK(cudaFree(d_out));
//     exit(-1);
// }
 
