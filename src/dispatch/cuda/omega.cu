#include <Rcpp.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>


#include "cutlass/numeric_conversion.h"

#include <cute/tensor.hpp>

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cuda_runtime.h>

#include <propr/context.h>
#include <propr/kernels/cuda/dispatch/omega.cuh>
#include <propr/kernels/cuda/detail/omega.cuh>


using namespace Rcpp;
using namespace propr;


void printRMatrix(Rcpp::NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      std::cout << mat(i, j) << " ";
    } std::cout << std::endl;
  } std::cout << std::endl;
}

NumericMatrix pad_matrix(const NumericMatrix& mat,
                         int padTop, int padBottom,
                         int padLeft, int padRight) {
  int old_nrow = mat.nrow();
  int old_ncol = mat.ncol();
  int new_nrow = old_nrow + padTop + padBottom;
  int new_ncol = old_ncol + padLeft + padRight;
  
  // create output matrix initialized to zero
  NumericMatrix out(new_nrow, new_ncol);
  
  // copy original into the right offset
  for (int i = 0; i < old_nrow; ++i) {
    for (int j = 0; j < old_ncol; ++j) {
      out(i + padTop, j + padLeft) = mat(i, j);
    }
  }
  
  return out;
}

void 
propr::dispatch::cuda::dof_global(NumericVector& out, const NumericMatrix& W, propr_context context) {
    using Config = kernels::cutlass_impl::OmegaConfig;
    NumericMatrix W2(W);
    // printRMatrix(W2);
    // for(int i=0; i < W2.nrow(); i++){
    //   for(int j=0; j < W2.ncol(); j++){
    //     W2(i,j) = float(335.712);
    //   } 
    // }

    for(int i=0; i < W2.nrow(); i++){
      for(int j=0; j < W2.ncol(); j++){
        W2(i,j) = float(0);
      } 
    }

    for(int i=0; i < W2.nrow(); i++){
      for(int j=0; j < W2.ncol(); j++){
        W2(i,j) = float(j) * W2.nrow() + i + 1;
      }
    }
    printRMatrix(W2);

    int t = 32;
    auto Wl = pad_matrix(W2, 0, ((W2.nrow() + t - 1)/t)*t - W2.nrow(), 0, ((W2.ncol() + t - 1)/t)*t - W2.ncol());
    int nfeats  = Wl.ncol();
    int samples = Wl.nrow();
    // printRMatrix(Wl);

    //  for(int i=0; i < Wl.nrow(); i++){
    //   for(int j=0; j < Wl.ncol(); j++){
    //     Wl(i,j) = float(i) + float(j)*Wl.nrow();
    //   }
    // }
    // printRMatrix(Wl);
    
    std::cout << "[" << W.nrow() << "," << W.ncol() << "]" << std::endl;

    const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    // CHECK_VECTOR_SIZE(out, llt);
    int alignment = 1; //16
    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP>(Wl, W_stride, alignment);
    float* d_out = nullptr;
    offset_t dout_stride = ((nfeats  + alignment - 1) / alignment) * alignment;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));

    const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    dim3 block(cute::size(Config::TiledMMA{}));
    dim3 grid(BX, BY);
    static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    const auto fptr = propr::kernels::cutlass_impl::omega_kernel<Config, cute::half_t, float>;

    cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, kShmSize);
    cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);

    std::cout <<  "M: "<< nfeats << " K: " << samples << std::endl;
    std::cout <<  "fptr<<<(" << grid.x << ","<< grid.y <<")"<< ",(" << block.x << "," << block.y <<")>>>" << std::endl;
    fptr<<<grid, block, kShmSize, context.stream>>>(nfeats, samples, d_W, d_out);
    
    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    auto h_full= new std::vector<float> (nfeats * dout_stride);
    CUDA_CHECK(cudaMemcpy(
        h_full->data(),
        d_out,
        nfeats * dout_stride * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    size_t counter = 0;
    double* out_ptr = REAL(out);
    std::cout << "[GPU]: \n";
    for (int i = 0; i < nfeats; ++i) {
        for (int j = 0; j < nfeats; ++j) {
            float v = h_full->at(size_t(i) * dout_stride + j);
            // out_ptr[counter++] = static_cast<double>(v);
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
    exit(-1);
}
 
void 
propr::dispatch::cuda::dof_population(NumericVector& out, const NumericMatrix& W, propr_context context) {
    // using Converter = cutlass::NumericConverter<cute::half_t, float>;
    // using Config     = kernels::cutlass_impl::OmegaConfig;

    // int nfeats   = W.ncol();
    // int samples = W.nrow();

    // const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    // CHECK_VECTOR_SIZE(out, llt);

    // offset_t W_stride;
    // auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP, true>(W, W_stride);

    // cute::half_t* d_out = nullptr;
    // CUDA_CHECK(cudaMalloc(&d_out, nfeats * nfeats * sizeof(cute::half_t)));

    // const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    // const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    // dim3 block(cute::size(Config::TiledMMA{}));
    // dim3 grid(BX, BY);
    // static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    // static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    // int shm_size = kShmSize;
    // const auto fptr = propr::kernels::cutlass_impl::omega_kernel<Config,cute::half_t>;

    // cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, shm_size);
    // cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
    // fptr<<<grid, block, shm_size, context.stream>>>(d_W, d_out, samples, W_stride);
    // CUDA_CHECK(cudaStreamSynchronize(context.stream));

    // std::vector<cute::half_t> h_full(nfeats * nfeats);
    // CUDA_CHECK(cudaMemcpy(
    //     h_full.data(),
    //     d_out,
    //     W_stride * W_stride * sizeof(cute::half_t),
    //     cudaMemcpyDeviceToHost
    // ));
    // Converter cast_away;
    // size_t counter = 0;
    // double* out_ptr = REAL(out);
    // for (int i = 1; i < nfeats; ++i) {
    //     for (int j = 0; j < i; ++j) {
    //         cute::half_t v = h_full[size_t(i) * W_stride + j];
    //         out_ptr[counter++] = static_cast<double>(cast_away(v));
    //     }
    // }

    // CUDA_CHECK(cudaFree(d_W));
    // CUDA_CHECK(cudaFree(d_out));
}