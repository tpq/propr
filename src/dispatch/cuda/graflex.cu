#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/transform_reduce.h>

#include <propr/utils.hpp>
#include <propr/kernels/cuda/dispatch/graflex.cuh>


using namespace Rcpp;
using namespace propr;

void
dispatch::cuda::getOR( NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getOR.");
    int a = 0, b = 0, c = 0, d = 0;
    double odds_ratio = static_cast<double>(a * d) / (b * c);
    double log_odds_ratio = std::log(odds_ratio);
}

void
dispatch::cuda::getORperm(NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getORperm.");
    int a = 0, b = 0, c = 0, d = 0;
    double odds_ratio = static_cast<double>(a * d) / (b * c);
    double log_odds_ratio = std::log(odds_ratio);
}

void
dispatch::cuda::permuteOR(NumericMatrix& out, const IntegerMatrix& A, const IntegerMatrix& G, int p, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for permuteOR.");
}

void
dispatch::cuda::getFDR(List& out, double actual, const NumericVector& permuted, propr::propr_context context) {
    const int n = permuted.size();
    double* d_permuted;
    CUDA_CHECK(cudaMalloc(&d_permuted, n * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(
        d_permuted, permuted.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_permuted(d_permuted);
    auto policy = thrust::cuda::par.on(context.stream);
    constexpr int2 init{0,0};
    int2 result = thrust::transform_reduce(
        policy,
        dev_permuted,  dev_permuted + n,
        [=] __host__ __device__ (double x) {
            return int2{ x >= actual, x <= actual };
        },
        init,
        [] __host__ __device__ (const int2& a, const int2& b) {
            return int2{ a.x + b.x, a.y + b.y };
        }
    );

    CUDA_CHECK(cudaStreamSynchronize(context.stream));
    CUDA_CHECK(cudaFree(d_permuted));

    double fdr_over = static_cast<double>(result.x) / n;
    double fdr_under = static_cast<double>(result.y) / n;
    out["over"] = fdr_over;
    out["under"] = fdr_under;
}

void
dispatch::cuda::getG(IntegerMatrix& out, const IntegerVector& Gk, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for getG.");
}

void
dispatch::cuda::graflex(NumericVector& out, const IntegerMatrix& A, const IntegerVector& Gk, int p, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for graflex.");   
}