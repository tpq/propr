#include <Rcpp.h>
#include <cuda_runtime.h>
#include <propr/interface/backend.hpp>

#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/rcpp_cuda.cuh>
#include <propr/utils/cuda_helpers.cuh>


#include <propr/kernels/cuda/dispatch/lrm.cuh>
#include <propr/kernels/cuda/detail/lrm.cuh>

using namespace Rcpp;

void
propr::dispatch::cuda::lrm_basic(NumericVector& out, NumericMatrix &Y, propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y;
    offset_t stride; d_Y = RcppMatrixToDevice<float>(Y, stride);

    float* d_mean;
    PROPR_CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim(propr::ceil_div(N_genes, blockDim.x),  propr::ceil_div(N_genes, blockDim.y));

    propr::detail::cuda::lrm_basic<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, stride, d_mean, N_samples, N_genes
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_mean, out, N_pairs);
    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_mean));
}

void
propr::dispatch::cuda::lrm_weighted(NumericVector& out,
                                    NumericMatrix &Y,
                                    NumericMatrix &W,
                                    propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

    
    offset_t stride_Y; float* d_Y = RcppMatrixToDevice<float>(Y, stride_Y);
    offset_t stride_W; float* d_W = RcppMatrixToDevice<float>(W, stride_W);


    float* d_mean;
    PROPR_CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim(propr::ceil_div(N_genes, blockDim.x), propr::ceil_div(N_genes, blockDim.y));

    propr::detail::cuda::lrm_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, stride_Y, d_W, stride_W, d_mean, N_samples, N_genes
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

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
    int N1      = Y.nrow();
    int N_genes = Y.ncol();
    int NT      = Yfull.nrow();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y; float* d_Yfull;
    offset_t stride_Y    ; d_Y     = RcppMatrixToDevice<float>(Y, stride_Y);
    offset_t stride_Yfull; d_Yfull = RcppMatrixToDevice<float>(Yfull, stride_Yfull);

    float* d_means;
    PROPR_CUDA_CHECK(cudaMalloc(&d_means, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim(propr::ceil_div(N_genes, blockDim.x), propr::ceil_div(N_genes, blockDim.y));

    propr::detail::cuda::lrm_alpha<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, stride_Y, d_Yfull, stride_Yfull, N1, NT, static_cast<float>(a), d_means, N1, N_genes
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_means, out, N_pairs);
    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_Yfull));
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

    dim3 blockDim(16, 16);
    dim3 gridDim(propr::ceil_div(N_genes, blockDim.x), propr::ceil_div(N_genes, blockDim.y));

    propr::detail::cuda::lrm_alpha_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y    , stride_Y, 
        d_Yfull, stride_Yfull,
        d_W    , stride_W,
        d_Wfull, stride_Wfull, 
        N1, NT, static_cast<float>(a), d_means, N_genes
    );
    PROPR_CUDA_CHECK(cudaGetLastError());
    PROPR_CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_means, out, N_pairs);

    PROPR_CUDA_CHECK(cudaFree(d_Y));
    PROPR_CUDA_CHECK(cudaFree(d_W));
    PROPR_CUDA_CHECK(cudaFree(d_Yfull));
    PROPR_CUDA_CHECK(cudaFree(d_Wfull));
    PROPR_CUDA_CHECK(cudaFree(d_means));
}