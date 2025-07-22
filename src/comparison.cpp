#include <Rcpp.h>
#include <propr/interface/comparsison.hpp>
#include <propr/interface/device_selector.hpp>

#include <propr/context.h>
#include <propr/kernels/cpu/dispatch/comparison.hpp>
#include <propr/kernels/cuda/dispatch/comparison.cuh>

using namespace propr;

// [[Rcpp::export]]
int count_less_than(Rcpp::NumericVector x, double cutoff, bool use_gpu) {
    if (is_gpu_backend() || use_gpu) {
        return dispatch::cuda::count_less_than(x, cutoff);
    } else {
        return dispatch::cpu::count_less_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_greater_than(Rcpp::NumericVector x, double cutoff, bool use_gpu) {
    if (is_gpu_backend() || use_gpu) {
        auto res = dispatch::cuda::count_greater_than(x, cutoff);
        return res;
    } else {
        return dispatch::cpu::count_greater_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_less_equal_than(Rcpp::NumericVector x, double cutoff, bool use_gpu) {
    if (is_gpu_backend() || use_gpu) {
        auto res = dispatch::cuda::count_less_equal_than(x, cutoff);
        return res;
    } else {
        return dispatch::cpu::count_less_equal_than(x, cutoff);
    }
}

// [[Rcpp::export]]
int count_greater_equal_than(Rcpp::NumericVector x, double cutoff, bool use_gpu) {
    if (is_gpu_backend() || use_gpu) {
        auto res = dispatch::cuda::count_greater_equal_than(x, cutoff);
        return res;
    } else {
        return dispatch::cpu::count_greater_equal_than(x, cutoff);
    }
}