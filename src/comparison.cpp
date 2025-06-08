#include <Rcpp.h>
#include <propr/interface/comparsison.hpp>
#include <propr/interface/device_selector.hpp>

#include <propr/context.h>
#include <propr/kernels/cpu/dispatch/comparison.hpp>
#include <propr/kernels/cuda/dispatch/comparison.cuh>

using namespace propr;

// [[Rcpp::export]]
int count_less_than(Rcpp::NumericVector x, double cutoff) {
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for count_less_than: %s. Falling back to CPU.", cudaGetErrorString(err));
            return dispatch::cpu::count_less_than(x, cutoff);
        } else {
            context.stream = stream;
            auto res = dispatch::cuda::count_less_than(x, cutoff,context);
            cudaStreamDestroy(stream);
            return res;
        }
    } else {
        return dispatch::cpu::count_less_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_greater_than(Rcpp::NumericVector x, double cutoff) {
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for count_greater_than: %s. Falling back to CPU.", cudaGetErrorString(err));
            return dispatch::cpu::count_greater_than(x, cutoff);
        } else {
            context.stream = stream;
            auto res = dispatch::cuda::count_greater_than(x, cutoff,context);
            cudaStreamDestroy(stream);
            return res;
        }
    } else {
        return dispatch::cpu::count_greater_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_less_equal_than(Rcpp::NumericVector x, double cutoff) {
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for count_less_equal_than: %s. Falling back to CPU.", cudaGetErrorString(err));
            return dispatch::cpu::count_less_equal_than(x, cutoff);
        } else {
            context.stream = stream;
            auto res = dispatch::cuda::count_less_equal_than(x, cutoff,context);
            cudaStreamDestroy(stream);
            return res;
        }
    } else {
        return dispatch::cpu::count_less_equal_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_greater_equal_than(Rcpp::NumericVector x, double cutoff) {
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for count_greater_equal_than: %s. Falling back to CPU.", cudaGetErrorString(err));
            return dispatch::cpu::count_greater_equal_than(x, cutoff);
        } else {
            context.stream = stream;
            auto res = dispatch::cuda::count_greater_equal_than(x, cutoff,context);
            cudaStreamDestroy(stream);
            return res;
        }
    } else {
        return dispatch::cpu::count_greater_equal_than(x, cutoff);
    }
}