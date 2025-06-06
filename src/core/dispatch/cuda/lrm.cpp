#include <Rcpp.h>
#include <cuda_runtime.h>
#include "../../../include/interface/backend.hpp"
#include "../../../include/utils.hpp"
#include "../../../include/kernels/cuda/dispatch/lrm.cuh"
#include "../../../include/kernels/cuda/detail/lrm.cuh"

using namespace Rcpp;

void
propr::dispatch::cuda::lrm_basic(NumericMatrix &Y, NumericVector& out, propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y;
    int num_rows, num_cols;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows, num_cols);

    float* d_mean;
    CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrm_basic<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_mean, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_mean, out, N_pairs);
    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_mean));
}

void
propr::dispatch::cuda::lrm_weighted(NumericMatrix &Y,
                                   NumericMatrix &W,
                                   NumericVector& out,
                                   propr::propr_context context) {
    int N_samples = Y.nrow();
    int N_genes = Y.ncol();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y;
    int num_rows_Y, num_cols_Y;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows_Y, num_cols_Y);

    float* d_W;
    int num_rows_W, num_cols_W;
    d_W = RcppNumericMatrixToDeviceFloat(W, num_rows_W, num_cols_W);


    float* d_mean;
    CUDA_CHECK(cudaMalloc(&d_mean, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrm_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_W, d_mean, N_samples, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_mean, out, N_pairs);
    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_mean));
}

void
propr::dispatch::cuda::lrm_alpha(NumericMatrix &Y,
                                     const double a,
                                     NumericMatrix& Yfull,
                                     NumericVector& out,
                                     propr::propr_context context) {
    int N1      = Y.nrow();
    int N_genes = Y.ncol();
    int NT      = Yfull.nrow();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);
    float* d_Y;
    int num_rows_Y, num_cols_Y;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows_Y, num_cols_Y);

    float* d_Yfull;
    int num_rows_Yfull, num_cols_Yfull;
    d_Yfull = RcppNumericMatrixToDeviceFloat(Yfull, num_rows_Yfull, num_cols_Yfull);

    float* d_means;
    CUDA_CHECK(cudaMalloc(&d_means, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrm_alpha<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_Yfull, N1, NT, static_cast<float>(a), d_means, N1, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_means, out, N_pairs);
    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_means));
}

void
propr::dispatch::cuda::lrm_alpha_weighted(NumericMatrix &Y,
                                              NumericMatrix &W,
                                              const double a,
                                              NumericMatrix& Yfull,
                                              NumericMatrix& Wfull,
                                              NumericVector& out,
                                              propr::propr_context context) {
    int N1 = Y.nrow();
    int N_genes = Y.ncol();
    int NT = Yfull.nrow();
    int N_pairs = N_genes * (N_genes - 1) / 2;
    CHECK_VECTOR_SIZE(out, N_pairs);

    float* d_Y;
    int num_rows_Y, num_cols_Y;
    d_Y = RcppNumericMatrixToDeviceFloat(Y, num_rows_Y, num_cols_Y);

    float* d_W;
    int num_rows_W, num_cols_W;
    d_W = RcppNumericMatrixToDeviceFloat(W, num_rows_W, num_cols_W);

    float* d_Yfull;
    int num_rows_Yfull, num_cols_Yfull;
    d_Yfull = RcppNumericMatrixToDeviceFloat(Yfull, num_rows_Yfull, num_cols_Yfull);

    float* d_Wfull;
    int num_rows_Wfull, num_cols_Wfull;
    d_Wfull = RcppNumericMatrixToDeviceFloat(Wfull, num_rows_Wfull, num_cols_Wfull);

    float* d_means;
    CUDA_CHECK(cudaMalloc(&d_means, N_pairs * sizeof(float)));

    dim3 blockDim(16, 16);
    dim3 gridDim((N_genes + blockDim.x - 1) / blockDim.x, (N_genes + blockDim.y - 1) / blockDim.y);

    propr::detail::cuda::lrm_alpha_weighted<<<gridDim, blockDim, 0, context.stream>>>(
        d_Y, d_Yfull, d_W, d_Wfull, N1, NT, static_cast<float>(a), d_means, N_genes
    );
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaStreamSynchronize(context.stream));

    copyFloatToNumericVector(d_means, out, N_pairs);

    CUDA_CHECK(cudaFree(d_Y));
    CUDA_CHECK(cudaFree(d_W));
    CUDA_CHECK(cudaFree(d_Yfull));
    CUDA_CHECK(cudaFree(d_Wfull));
    CUDA_CHECK(cudaFree(d_means));
}