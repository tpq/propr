#include <Rcpp.h>

#include <propr/interface/comparsison.hpp>
#include <propr/kernels/cpu/dispatch/comparison.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
int count_less_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto") {
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        return runtime::cuda_executor::count_less_than(x, cutoff);
    }
    return dispatch::cpu::count_less_than(x, cutoff);
}

// [[Rcpp::export]]
int count_greater_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto") {
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        return runtime::cuda_executor::count_greater_than(x, cutoff);
    }
    return dispatch::cpu::count_greater_than(x, cutoff);
}

// [[Rcpp::export]]
int count_less_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto") {
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        return runtime::cuda_executor::count_less_equal_than(x, cutoff);
    }
    return dispatch::cpu::count_less_equal_than(x, cutoff);
}

// [[Rcpp::export]]
int count_greater_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto") {
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        return runtime::cuda_executor::count_greater_equal_than(x, cutoff);
    }
    return dispatch::cpu::count_greater_equal_than(x, cutoff);
}
