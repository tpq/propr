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

struct threshold_counter_handle {
    runtime::Backend backend;
    void* impl;

    ~threshold_counter_handle() {destroy(); }

    void destroy() {
        if (impl == nullptr) return;
        if (backend == runtime::Backend::CUDA) {
            runtime::cuda_executor::count_values_beyond_thresholds_destroy(impl);
        } else {
            dispatch::cpu::count_values_beyond_thresholds_destroy(impl);
        }
        impl = nullptr;
    }
};

// [[Rcpp::export]]
SEXP count_values_beyond_thresholds_begin(
    Rcpp::NumericVector cutoffs,
    bool direct,
    Rcpp::String backend = "auto") {
    runtime::Backend resolved = runtime::resolve_backend(backend);
    auto* handle = new threshold_counter_handle{resolved, nullptr};

    if (resolved == runtime::Backend::CUDA) {
        handle->impl =runtime::cuda_executor::count_values_beyond_thresholds_begin( cutoffs, direct);
    } else {
        handle->impl = dispatch::cpu::count_values_beyond_thresholds_begin(cutoffs, direct);
    }

    Rcpp::XPtr<threshold_counter_handle> ptr(handle, true);
    return ptr;
}

// [[Rcpp::export]]
void count_values_beyond_thresholds_accumulate(
    SEXP counter,
    Rcpp::NumericVector values) {
    Rcpp::XPtr<threshold_counter_handle> ptr(counter);

    if (ptr->backend == runtime::Backend::CUDA) {
        runtime::cuda_executor::count_values_beyond_thresholds_accumulate( ptr->impl, values);
    } else {
        dispatch::cpu::count_values_beyond_thresholds_accumulate(ptr->impl, values);
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector count_values_beyond_thresholds_end(SEXP counter) {
    Rcpp::XPtr<threshold_counter_handle> ptr(counter);
    threshold_counter_handle* raw = ptr.get();

    Rcpp::NumericVector out;
    if (raw->backend == runtime::Backend::CUDA) {
        out = runtime::cuda_executor::count_values_beyond_thresholds_end(raw->impl);
    } else {
        out = dispatch::cpu::count_values_beyond_thresholds_end(raw->impl);
    }

    raw->destroy();
    delete raw;
    R_ClearExternalPtr(counter);
    return out;
}
