#include <cfloat>
#include <math.h>
#include <stdio.h>
#include "../../../include/utils.hpp"
#include "../../../include/kernels/cuda/detail/lrv.cuh"
#include "../../../include/kernels/cuda/dispatch/lrv.cuh"
using namespace Rcpp;



void
propr::dispatch::cuda::lrv_basic(Rcpp::NumericMatrix &Y,
                                 Rcpp::NumericVector& out,
                                 propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs); // Check output vector size

    float* d_Y;
    int num_rows, num_cols;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows, num_cols);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrv_basic<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_variances, out, N_pairs);
    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_variances));
}

void
propr::dispatch::cuda::lrv_weighted(Rcpp::NumericMatrix &Y,
                                    Rcpp::NumericMatrix &W,
                                    Rcpp::NumericVector& out,
                                    propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs); // Check output vector size

    float* d_Y;
    float* d_W;
    int num_rows, num_cols;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows, num_cols);
    d_W = RcppNumericMatrixToDeviceFloat(W, num_rows, num_cols);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrv_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_W, d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_variances, out, N_pairs); // Copy result to 'out' parameter

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_variances));
}

void propr::dispatch::cuda::lrv_alpha(Rcpp::NumericMatrix &Y,
                                           double a,
                                           Rcpp::NumericMatrix& Yfull,
                                           Rcpp::NumericVector& out,
                                           propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y;
    float* d_Yfull;
    int num_rows, num_cols;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows, num_cols);
    d_Yfull = RcppNumericMatrixToDeviceFloat(Yfull, num_rows, num_cols);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrv_alpha<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_Yfull, static_cast<float>(a), d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_variances, out, N_pairs);

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_variances));
}

void
propr::dispatch::cuda::lrv_alpha_weighted(Rcpp::NumericMatrix &Y,
                                          Rcpp::NumericMatrix &W,
                                          double a,
                                          Rcpp::NumericMatrix& Yfull,
                                          Rcpp::NumericMatrix& Wfull,
                                          Rcpp::NumericVector& out,
                                          propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y,* d_W,* d_Yfull,* d_Wfull;
    int num_rows, num_cols;
    d_Y     = RcppNumericMatrixToDeviceFloat(Y, num_rows, num_cols);
    d_W     = RcppNumericMatrixToDeviceFloat(W, num_rows, num_cols);
    d_Yfull = RcppNumericMatrixToDeviceFloat(Yfull, num_rows, num_cols);
    d_Wfull = RcppNumericMatrixToDeviceFloat(Wfull, num_rows, num_cols);

    float* d_variances;
    CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrv_alpha_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_Yfull, d_W, d_Wfull, static_cast<float>(a), d_variances, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_variances, out, N_pairs); // Copy result to 'out' parameter

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_Wfull));
    CUDA_CHECK(cudaFree(d_variances));
}