#include <Rcpp.h>

#include <propr/interface/omega.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/omega.hpp>
#include <propr/kernels/cuda/dispatch/omega.cuh>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector omega(Rcpp::NumericMatrix & W, bool use_gpu) {
    size_t nfeats = W.ncol();
    size_t llt    = nfeats * (nfeats - 1 ) / 2;
    Rcpp::NumericVector result(llt);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::dof_global(result, W);
    } else {
        dispatch::cpu::dof_global(result, W);
    }
    return result;
}


// [[Rcpp::export]]
Rcpp::NumericVector Omega(Rcpp::NumericMatrix & W, bool use_gpu) {
    size_t nfeats = W.ncol();
    size_t llt    = nfeats * (nfeats - 1 ) / 2;
    Rcpp::NumericVector result(llt);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::dof_population(result, W);
    } else {
        dispatch::cpu::dof_population(result, W);
    }
    return result;
}