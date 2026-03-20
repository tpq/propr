#include <Rcpp.h>
#include <cuda_runtime.h>
#include <propr/interface/backend.hpp>

#include <propr/utils/rcpp/rcpp_checks.h>
#include <propr/utils/cuda/cuda_checks.h>
#include <propr/utils/rcpp/rcpp_cuda.cuh>
#include <propr/utils/common/cuda_helpers.cuh>
#include <propr/utils/profilers/cuda_profiler.cuh>


#include <propr/kernels/cuda/dispatch/lrm.cuh>
#include <propr/kernels/cuda/detail/lrm.cuh>
#include <propr/kernels/cuda/traits/lrm.cuh>


using namespace Rcpp;

void
propr::dispatch::cuda::lrm_basic(NumericVector& out, NumericMatrix &Y, propr::propr_context context) {
    using Config = propr::cuda::traits::lrm_basic;
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y;
    offset_t stride; d_Y = RcppMatrixToDevice<float>(Y, stride);

    float* d_mean_log;
    PROPR_CUDA_CHECK(cudaMalloc(&d_mean_log, static_cast<size_t>(N_genes) * sizeof(float)));

    float* d_mean;
    PROPR_CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 block1(Config::P1_Layout::BLK_X);
    dim3 grid1(propr::ceil_div(N_genes, Config::P1_Layout::BLK_X));

    dim3 block2(Config::P2_Layout::BLK_X, Config::P2_Layout::BLK_Y);
    dim3 grid2(propr::ceil_div(N_genes, Config::P2_Layout::BLK_X),
               propr::ceil_div(N_genes, Config::P2_Layout::BLK_Y));
               

    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        propr::detail::cuda::lrm_basic_phase_1<Config><<<grid1, block1, 0, context.stream>>>(
            d_Y, stride, d_mean_log, N_samples, N_genes
        );
        PROPR_CUDA_CHECK(cudaGetLastError());

        propr::detail::cuda::lrm_basic_phase_2<Config><<<grid2, block2, 0, context.stream>>>(
            d_mean_log, d_mean, N_genes
        );
        PROPR_CUDA_CHECK(cudaGetLastError());
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    copyToNumericVector(d_mean, out, N_pairs);
    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_mean_log));
    PROPR_CUDA_CHECK(cudaFree(d_mean));
}

void
propr::dispatch::cuda::lrm_weighted(NumericVector& out,
                                    NumericMatrix &Y,
                                    NumericMatrix &W,
                                    propr::propr_context context) {
    using Config = propr::cuda::traits::lrm_weighted;
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

    
    offset_t stride_Y; float* d_Y = RcppMatrixToDevice<float>(Y, stride_Y);
    offset_t stride_W; float* d_W = RcppMatrixToDevice<float>(W, stride_W);


    float* d_mean;
    PROPR_CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 blockDim(Config::BLK_X, Config::BLK_Y);
    dim3 gridDim(propr::ceil_div(N_genes, Config::BLK_X), propr::ceil_div(N_genes, Config::BLK_Y));

    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        propr::detail::cuda::lrm_weighted<Config><<<gridDim, blockDim, 0, context.stream>>>(
            d_Y, stride_Y, d_W, stride_W, d_mean, N_samples, N_genes
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    copyToNumericVector(d_mean, out, N_pairs);
    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_W));
    PROPR_CUDA_CHECK(cudaFree(d_mean));
}

void
propr::dispatch::cuda::lrm_alpha(NumericVector& out,
                                 NumericMatrix &Y,
                                 const double a,
                                 NumericMatrix& Yfull,
                                 propr::propr_context context) {
    using Config = propr::cuda::traits::lrm_alpha;

    int N1      = Y.nrow();
    int N_genes = Y.ncol();
    int NT      = Yfull.nrow();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y; float* d_Yfull;
    offset_t stride_Y    ; d_Y     = RcppMatrixToDevice<float>(Y, stride_Y);
    offset_t stride_Yfull; d_Yfull = RcppMatrixToDevice<float>(Yfull, stride_Yfull);

    float* d_h;
    PROPR_CUDA_CHECK(cudaMalloc(&d_h, static_cast<size_t>(N_genes) * sizeof(float)));

    float* d_means;
    PROPR_CUDA_CHECK(cudaMalloc(&d_means, N_pairs * sizeof(float)));

    dim3 block1(Config::P1_Layout::BLK_X);
    dim3 grid1(propr::ceil_div(N_genes, Config::P1_Layout::BLK_X));

    dim3 block2(Config::P2_Layout::BLK_X, Config::P2_Layout::BLK_Y);
    dim3 grid2(propr::ceil_div(N_genes, Config::P2_Layout::BLK_X),
               propr::ceil_div(N_genes, Config::P2_Layout::BLK_Y));

    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        propr::detail::cuda::lrm_alpha_phase_1<Config><<<grid1, block1, 0, context.stream>>>(
            d_Y, stride_Y, d_Yfull, stride_Yfull, N1, NT, static_cast<float>(a), d_h, N_genes
        );
        PROPR_CUDA_CHECK(cudaGetLastError());

        propr::detail::cuda::lrm_alpha_phase_2<Config><<<grid2, block2, 0, context.stream>>>(
            d_h, d_means, N_genes
        );
        PROPR_CUDA_CHECK(cudaGetLastError());
        PROPR_STREAM_SYNCHRONIZE(context);

    }
    copyToNumericVector(d_means, out, N_pairs);
    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_Yfull));
    PROPR_CUDA_CHECK(cudaFree(d_h));
    PROPR_CUDA_CHECK(cudaFree(d_means));
}

void
propr::dispatch::cuda::lrm_alpha_weighted(NumericVector& out,
                                          NumericMatrix &Y,
                                          NumericMatrix &W,
                                          const double a,
                                          NumericMatrix& Yfull,
                                          NumericMatrix& Wfull,
                                          propr::propr_context context) {
    using Config = propr::cuda::traits::lrm_alpha_weighted;

    int N1 = Y.nrow();
    int N_genes = Y.ncol();
    int NT = Yfull.nrow();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y, * d_W, * d_Yfull, * d_Wfull;
    offset_t stride_Y    ; d_Y     = RcppMatrixToDevice<float>(Y, stride_Y); 
    offset_t stride_W    ; d_W     = RcppMatrixToDevice<float>(W, stride_W);
    offset_t stride_Yfull; d_Yfull = RcppMatrixToDevice<float>(Yfull, stride_Yfull);
    offset_t stride_Wfull; d_Wfull = RcppMatrixToDevice<float>(Wfull, stride_Wfull);

    float* d_means;
    PROPR_CUDA_CHECK(cudaMalloc(&d_means, N_pairs * sizeof(float)));

    dim3 blockDim(Config::BLK_X, Config::BLK_Y);
    dim3 gridDim(propr::ceil_div(N_genes, Config::BLK_X), propr::ceil_div(N_genes, Config::BLK_Y));

    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        propr::detail::cuda::lrm_alpha_weighted<Config><<<gridDim, blockDim, 0, context.stream>>>(
            d_Y    , stride_Y, 
            d_Yfull, stride_Yfull,
            d_W    , stride_W,
            d_Wfull, stride_Wfull, 
            N1, NT, static_cast<float>(a), d_means, N_genes
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    copyToNumericVector(d_means, out, N_pairs);

    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_W));
    PROPR_CUDA_CHECK(cudaFree(d_Yfull));
    PROPR_CUDA_CHECK(cudaFree(d_Wfull));
    PROPR_CUDA_CHECK(cudaFree(d_means));
}
