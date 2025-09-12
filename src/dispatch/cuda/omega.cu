#include <Rcpp.h>

#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cuda_runtime.h>

#include <propr/context.h>
#include <propr/kernels/cuda/dispatch/omega.cuh>
#include <propr/kernels/cuda/detail/omega.cuh>


#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>
#include <propr/utils/rcpp_helpers.h>

using namespace Rcpp;
using namespace propr;

void 
propr::dispatch::cuda::dof_global(NumericVector& out, const NumericMatrix& W, propr_context context) {
  using Config = omega_global_config;

    int t = 128;
    auto Wl = rcpp::helpers::pad_matrix(W, 0, ((W.nrow() + t - 1)/t)*t - W.nrow(), 0, ((W.ncol() + t - 1)/t)*t - W.ncol());
    int nfeats  = Wl.ncol();
    int samples = Wl.nrow();

    const size_t llt = static_cast<size_t>(W.ncol()) * (W.ncol() - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<float>(Wl, W_stride);
    
    float* d_out = nullptr;
    offset_t dout_stride = nfeats;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
    
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid(nfeats / Config::BLK_M, nfeats / Config::BLK_M);
    
    omega_kernel<Config><<<grid, block, 0, context.stream>>>(d_W, d_out, nfeats, samples);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full= new std::vector<float> (nfeats * nfeats);
    CUDA_CHECK(cudaMemcpy(
        h_full->data(),
        d_out,
        nfeats * nfeats * sizeof(float),
        cudaMemcpyDeviceToHost
    ));
    size_t counter = 0;
    double* out_ptr = REAL(out);
    for (int i = 1; i < W.ncol(); i++) {
        for (int j = 0; j < i; j++) {
            float v = h_full->at(size_t(i) * nfeats + j);
            out_ptr[counter++] = static_cast<double>(v);
        }
    }
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
}
 
void 
propr::dispatch::cuda::dof_population(NumericVector& out, const NumericMatrix& W, propr_context context) {
    using Config = omega_population_config;

    int t = 128;
    auto Wl = rcpp::helpers::pad_matrix(W, 0, ((W.nrow() + t - 1)/t)*t - W.nrow(), 0, ((W.ncol() + t - 1)/t)*t - W.ncol());
    int nfeats  = Wl.ncol();
    int samples = Wl.nrow();

    const size_t llt = static_cast<size_t>(W.ncol()) * (W.ncol() - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<float, REALSXP>(Wl, W_stride);
    
    float* d_out = nullptr; offset_t dout_stride = nfeats;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * dout_stride * sizeof(*d_out)));
    
    dim3 block(Config::BLK_M / Config::TH_X, Config::BLK_M / Config::TH_Y);
    dim3 grid(nfeats / Config::BLK_M, nfeats / Config::BLK_M);

    // std::cout << "fpt<<<(" << grid.x << "," << grid.y << ")" <<",(" << block.x << "," << block.y << ")>>>" << std::endl;
    omega_kernel<Config><<<grid, block, 0, context.stream>>>(d_W, d_out, nfeats, samples);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    auto h_full= new std::vector<float> (nfeats * nfeats);
    CUDA_CHECK(cudaMemcpy(
        h_full->data(),
        d_out,
        nfeats * nfeats * sizeof(float),
        cudaMemcpyDeviceToHost
    ));

    size_t counter = 0;
    double* out_ptr = REAL(out);
    for (int i = 1; i < W.ncol(); i++) {
        for (int j = 0; j < i; j++) {
            float v = h_full->at(size_t(i) * nfeats + j);
            out_ptr[counter++] = static_cast<double>(v);
        }
    }
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
}