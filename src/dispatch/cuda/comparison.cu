#include <Rcpp.h>

#include <thrust/count.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>

#include <propr/context.h>

#include <propr/kernels/cuda/dispatch/comparison.cuh>

#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>
#include <propr/utils/cuda_profiler.cuh>


using namespace Rcpp;
using namespace propr;


int 
dispatch::cuda::count_less_than(Rcpp::NumericVector& x,
                                double cutoff,
                                propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v < cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int 
dispatch::cuda::count_greater_than(Rcpp::NumericVector& x,
                                   double cutoff,
                                   propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v > cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }
    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int dispatch::cuda::count_less_equal_than(Rcpp::NumericVector& x,
                                          double cutoff,
                                          propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt =  0;
    {    
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v <= cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    }

    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}

int dispatch::cuda::count_greater_equal_than(Rcpp::NumericVector& x,
                                             double cutoff,
                                             propr::propr_context context) {
    const int n = x.size();
    double* d_x = nullptr;
    PROPR_CUDA_CHECK(cudaMalloc(&d_x, n * sizeof(double)));
    PROPR_CUDA_CHECK(cudaMemcpy(
        d_x, x.begin(), n * sizeof(double),
        cudaMemcpyHostToDevice
    ));

    thrust::device_ptr<double> dev_ptr(d_x);
    auto policy = thrust::cuda::par.on(context.stream);
    int cnt = 0;
    {
        PROPR_PROFILE_CUDA("kernel", context.stream);
        cnt = thrust::count_if(
            policy,
            dev_ptr, dev_ptr + n,
            [cutoff] __device__ (double v) { return v >= cutoff; }
        );
        PROPR_STREAM_SYNCHRONIZE(context);
    
    }
    PROPR_CUDA_CHECK(cudaFree(d_x));
    return cnt;
}
