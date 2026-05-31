#include <cfloat>
#include <math.h>
#include <stdio.h>

#include <propr/context.h>
#include <propr/runtime/dispatch.hpp>
#include <propr/utils/common/cuda_helpers.cuh>
#include <propr/utils/rcpp/rcpp_checks.h>
#include <propr/utils/rcpp/rcpp_cuda.cuh>
#include <propr/utils/profilers/cuda_profiler.cuh>


#include <propr/kernels/cuda/detail/lrv.cuh>
#include <propr/kernels/cuda/dispatch/lrv.cuh>
#include <propr/kernels/cuda/traits/lrv.cuh>


using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::lrv_basic(Rcpp::NumericVector& out,
                          Rcpp::NumericMatrix &Y,
                          propr_context context) {
    propr::runtime::with_precision([&](auto tag) {
        using Real = typename decltype(tag)::type;
        using Config = propr::cuda::traits::lrv_basic;
        int N_samples = Y.nrow();
        int N_genes   = Y.ncol();
        size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
        PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

        Real* d_Y;
        offset_t stride; d_Y = RcppMatrixToDevice<Real>(Y, stride);
        Real* d_variances;
        PROPR_CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(Real)));

        dim3 blockDim(Config::BLK_X, Config::BLK_Y);
        dim3 gridDim(propr::ceil_div(N_genes,Config::BLK_X), propr::ceil_div(N_genes,Config::BLK_Y));

        {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            detail::cuda::lrv_basic<Real, Config><<<gridDim, blockDim, 0, context.stream>>>(
                d_Y, stride, d_variances, N_samples, N_genes
            );
            PROPR_STREAM_SYNCHRONIZE(context);
        }

        copyToNumericVector(d_variances, out, N_pairs);
        PROPR_CUDA_CHECK(cudaFree(d_Y));
        PROPR_CUDA_CHECK(cudaFree(d_variances));
    });
}

void
dispatch::cuda::lrv_weighted(Rcpp::NumericVector& out,
                             Rcpp::NumericMatrix &Y,
                             Rcpp::NumericMatrix &W,
                             propr_context context) {
    propr::runtime::with_precision([&](auto tag) {
        using Real = typename decltype(tag)::type;
        using Config = propr::cuda::traits::lrv_weighted;
        int N_samples = Y.nrow();
        int N_genes = Y.ncol();
        size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
        PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

        Real* d_Y, * d_W;
        offset_t Y_stride; d_Y = RcppMatrixToDevice<Real>(Y, Y_stride);
        offset_t W_stride; d_W = RcppMatrixToDevice<Real>(W, W_stride);

        Real* d_variances;
        PROPR_CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(Real)));

        dim3 blockDim(Config::BLK_X, Config::BLK_Y);
        dim3 gridDim(propr::ceil_div(N_genes,Config::BLK_X), propr::ceil_div(N_genes,Config::BLK_Y));

       {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            detail::cuda::lrv_weighted<Real, Config><<<gridDim, blockDim, 0, context.stream>>>(
                d_Y, Y_stride,
                d_W, W_stride,
                d_variances, N_samples, N_genes
            );
            PROPR_STREAM_SYNCHRONIZE(context);
       }
        copyToNumericVector(d_variances, out, N_pairs);

        PROPR_CUDA_CHECK(cudaFree(d_Y));
        PROPR_CUDA_CHECK(cudaFree(d_W));
        PROPR_CUDA_CHECK(cudaFree(d_variances));
    });
}

void 
dispatch::cuda::lrv_alpha(Rcpp::NumericVector& out,
                          Rcpp::NumericMatrix &Y,
                          double a,
                          Rcpp::NumericMatrix& Yfull,
                          propr_context context) {
    propr::runtime::with_precision([&](auto tag) {
        using Real = typename decltype(tag)::type;
        using Config = propr::cuda::traits::lrv_alpha;
        int N_samples = Y.nrow();
        int N_samples_full = Yfull.nrow();
        int N_genes = Y.ncol();

        size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
        PROPR_CHECK_VECTOR_SIZE(out, N_pairs);

        Real* d_Y;
        Real* d_Yfull;

        offset_t Y_stride    ; d_Y     = RcppMatrixToDevice<Real>(Y,         Y_stride);
        offset_t Yfull_stride; d_Yfull = RcppMatrixToDevice<Real>(Yfull, Yfull_stride);

        Real* d_variances;
        PROPR_CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(Real)));

        dim3 blockDim(Config::BLK_X, Config::BLK_Y);
        dim3 gridDim(propr::ceil_div(N_genes, Config::BLK_X), propr::ceil_div(N_genes, Config::BLK_Y));

        {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            detail::cuda::lrv_alpha<Real, Config><<<gridDim, blockDim, 0, context.stream>>>(
                d_Y    , Y_stride,
                d_Yfull, Yfull_stride,
                static_cast<Real>(a), d_variances, N_samples, N_samples_full, N_genes
            );
        }

        PROPR_STREAM_SYNCHRONIZE(context);

        copyToNumericVector(d_variances, out, N_pairs);

        PROPR_CUDA_CHECK(cudaFree(d_Y));
        PROPR_CUDA_CHECK(cudaFree(d_Yfull));
        PROPR_CUDA_CHECK(cudaFree(d_variances));
    });
}

void
dispatch::cuda::lrv_alpha_weighted(Rcpp::NumericVector& out,
                                   Rcpp::NumericMatrix &Y,
                                   Rcpp::NumericMatrix &W,
                                   double a,
                                   Rcpp::NumericMatrix& Yfull,
                                   Rcpp::NumericMatrix& Wfull,
                                   propr_context context) {
    propr::runtime::with_precision([&](auto tag) {
        using Real = typename decltype(tag)::type;
        using Config = propr::cuda::traits::lrv_alpha_weighted;

        int N_samples = Y.nrow();
        int N_samples_full = Yfull.nrow();
        int N_genes = Y.ncol();
        size_t N_pairs = size_t(N_genes) * (N_genes - 1) / 2;
        PROPR_CHECK_VECTOR_SIZE(out, N_pairs);
        Real* d_Y,* d_W,* d_Yfull,* d_Wfull;
        offset_t Y_stride    ; d_Y     = RcppMatrixToDevice<Real>(Y    , Y_stride    );
        offset_t W_stride    ; d_W     = RcppMatrixToDevice<Real>(W    , W_stride    );
        offset_t Yfull_stride; d_Yfull = RcppMatrixToDevice<Real>(Yfull, Yfull_stride);
        offset_t Wfull_stride; d_Wfull = RcppMatrixToDevice<Real>(Wfull, Wfull_stride);

        Real* d_variances;
        PROPR_CUDA_CHECK(cudaMalloc(&d_variances, N_pairs * sizeof(Real)));

        dim3 blockDim(Config::BLK_X, Config::BLK_Y);
        dim3 gridDim(propr::ceil_div(N_genes,Config::BLK_X), propr::ceil_div(N_genes,Config::BLK_Y));

        {
            PROPR_PROFILE_CUDA("kernel", context.stream);
            detail::cuda::lrv_alpha_weighted<Real, Config><<<gridDim, blockDim, 0, context.stream>>>(
                d_Y    , Y_stride,
                d_Yfull, Yfull_stride,
                d_W    , W_stride,
                d_Wfull, Wfull_stride,
                static_cast<Real>(a), d_variances, N_samples, N_samples_full, N_genes
            );
            PROPR_STREAM_SYNCHRONIZE(context);
        }

        copyToNumericVector(d_variances, out, N_pairs);

        PROPR_CUDA_CHECK(cudaFree(d_Y));
        PROPR_CUDA_CHECK(cudaFree(d_W));
        PROPR_CUDA_CHECK(cudaFree(d_Yfull));
        PROPR_CUDA_CHECK(cudaFree(d_Wfull));
        PROPR_CUDA_CHECK(cudaFree(d_variances));
    });
}
