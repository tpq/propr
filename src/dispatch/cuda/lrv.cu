#include <cfloat>
#include <math.h>
#include <stdio.h>

#include <propr/context.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/rcpp_cuda.cuh>


#include <propr/kernels/cuda/detail/lrv.cuh>
#include <propr/kernels/cuda/dispatch/lrv.cuh>


using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::lrv_basic(Rcpp::NumericVector& out,
                          Rcpp::NumericMatrix &Y,
                          propr_context context) {
    int N_samples = Y.nrow();
    int N_genes   = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y;
    offset_t stride; d_Y = RcppMatrixToDevice<float>(Y, stride);
    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    detail::cuda::lrv_basic<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, stride, d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_variances, out, N_pairs);
    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_variances));
}

void
dispatch::cuda::lrv_weighted(Rcpp::NumericVector& out,
                             Rcpp::NumericMatrix &Y,
                             Rcpp::NumericMatrix &W,
                             propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y, * d_W;
    offset_t Y_stride; d_Y = RcppMatrixToDevice<float>(Y, Y_stride);
    offset_t W_stride; d_W = RcppMatrixToDevice<float>(W, W_stride);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    detail::cuda::lrv_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, Y_stride, 
        d_W, W_stride, 
        d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_variances, out, N_pairs);

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_variances));
}

void dispatch::cuda::lrv_alpha(Rcpp::NumericVector& out,
                               Rcpp::NumericMatrix &Y,
                               double a,
                               Rcpp::NumericMatrix& Yfull,
                               propr_context context) {
    int N_samples = Y.nrow();
    int N_samples_full = Yfull.nrow();
    int N_genes = Y.ncol();

    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y;
    float* d_Yfull;

    offset_t Y_stride    ; d_Y     = RcppMatrixToDevice<float>(Y,         Y_stride);
    offset_t Yfull_stride; d_Yfull = RcppMatrixToDevice<float>(Yfull, Yfull_stride);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    detail::cuda::lrv_alpha<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y    , Y_stride,
        d_Yfull, Yfull_stride,
        static_cast<float>(a), d_variances, N_samples, N_samples_full, N_genes
    );
    
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_variances, out, N_pairs);

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_variances));
}

void
dispatch::cuda::lrv_alpha_weighted(Rcpp::NumericVector& out,
                                   Rcpp::NumericMatrix &Y,
                                   Rcpp::NumericMatrix &W,
                                   double a,
                                   Rcpp::NumericMatrix& Yfull,
                                   Rcpp::NumericMatrix& Wfull,
                                   propr_context context) {
    int N_samples = Y.nrow();
    int N_samples_full = Yfull.nrow();
    int N_genes = Y.ncol();
    size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y,* d_W,* d_Yfull,* d_Wfull;
    offset_t Y_stride    ; d_Y     = RcppMatrixToDevice<float>(Y    , Y_stride    );
    offset_t W_stride    ; d_W     = RcppMatrixToDevice<float>(W    , W_stride    );
    offset_t Yfull_stride; d_Yfull = RcppMatrixToDevice<float>(Yfull, Yfull_stride);
    offset_t Wfull_stride; d_Wfull = RcppMatrixToDevice<float>(Wfull, Wfull_stride);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    detail::cuda::lrv_alpha_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y    , Y_stride, 
        d_Yfull, Yfull_stride, 
        d_W    , W_stride,
        d_Wfull, Wfull_stride,
        static_cast<float>(a), d_variances, N_samples, N_samples_full, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyToNumericVector(d_variances, out, N_pairs);

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_Wfull));
    CUDA_CHECK(cudaFree(d_variances));
}