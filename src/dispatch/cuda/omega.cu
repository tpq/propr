#include <Rcpp.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>

#include <propr/context.h>
#include <propr/kernels/cuda/dispatch/omega.cuh>
#include <propr/kernels/cuda/detail/omega.cuh>

#include <propr/utils/rcpp_cuda.cuh>


using namespace Rcpp;
using namespace propr;


void dof_global(NumericVector& out,
                const NumericMatrix& W,
                propr_context context) {
    using Converter = cutlass::NumericConverter<cute::half_t, float>;
    using Config     = kernels::cutlass::OmegaConfig;

    int nfeats = W.ncol();
    const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);

    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP, true>(W, W_stride);

    cute::half_t* d_out = nullptr;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats * nfeats * sizeof(cute::half_t)));

    const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    dim3 block(cute::size(Config::TiledMMA{}));
    dim3 grid(BX, BY);
    static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    int shm_size = kShmSize;
    const auto fptr = propr::kernels::cutlass::omega_kernel<Config,cute::half_t>;

    cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, shm_size);
    cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
    fptr<<<grid, block, shm_size, context.stream>>>(d_W, d_out, nfeats, nfeats);
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    std::vector<cute::half_t> h_full(nfeats * nfeats);
    CUDA_CHECK(cudaMemcpy(
        h_full.data(),
        d_out,
        nfeats * nfeats * sizeof(cute::half_t),
        cudaMemcpyDeviceToHost
    ));
    Converter cast_away;
    size_t counter = 0;
    double* out_ptr = REAL(out);
    for (int i = 1; i < nfeats; ++i) {
        for (int j = 0; j < i; ++j) {
            cute::half_t v = h_full[size_t(i) * nfeats + j];
            out_ptr[counter++] = static_cast<double>(cast_away(v));
        }
    }

    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
}

void dof_population(NumericVector& out,
                    const NumericMatrix& W,
                    propr_context context) {
    using Converter = cutlass::NumericConverter<cute::half_t, float>;
    using Config = kernels::cutlass::OmegaConfig;

    int nfeats = W.ncol();
    const size_t llt = static_cast<size_t>(nfeats) * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);

    offset_t W_stride;
    auto *d_W = RcppMatrixToDevice<cute::half_t, REALSXP, true>(W, W_stride);

    cute::half_t* d_out = nullptr;
    CUDA_CHECK(cudaMalloc(&d_out, nfeats* nfeats * sizeof(cute::half_t)));

    const int BX = (nfeats + Config::BLK_M - 1) / Config::BLK_M;
    const int BY = (nfeats + Config::BLK_M - 1) / Config::BLK_M;

    dim3 block(cute::size(Config::TiledMMA{}));
    dim3 grid(BX, BY);
    static constexpr int shm_size_AB = cute::cosize(Config::SmemLayoutA{}) + cute::cosize(Config::SmemLayoutB{});
    static constexpr int kShmSize = shm_size_AB * sizeof(__half);
    int shm_size = kShmSize;
    const auto fptr = propr::kernels::cutlass::omega_kernel<Config,cute::half_t>;

    cudaFuncSetAttribute(fptr, cudaFuncAttributeMaxDynamicSharedMemorySize, shm_size);
    cudaFuncSetAttribute(fptr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
    fptr<<<grid, block, shm_size, context.stream>>>(d_W, d_out, nfeats, nfeats);

    std::vector<cute::half_t> h_full(nfeats * nfeats);
    CUDA_CHECK(cudaMemcpy(
        h_full.data(),
        d_out,
        nfeats * nfeats * sizeof(cute::half_t),
        cudaMemcpyDeviceToHost
    ));
    Converter cast_away;
    size_t counter = 0;
    double* out_ptr = REAL(out);
    for (int i = 1; i < nfeats; ++i) {
        for (int j = 0; j < i; ++j) {
            cute::half_t v = h_full[size_t(i) * nfeats + j];
            out_ptr[counter++] = static_cast<double>(cast_away(v));
        }
    }


    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_out));
}