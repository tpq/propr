#include <Rcpp.h>
#include <propr/interface/ctzRcpp.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/ctzRcpp.hpp>
#include <propr/kernels/cuda/dispatch/ctzRcpp.cuh>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector ctzRcpp(Rcpp::NumericMatrix & X, bool use_gpu) {
    int nfeats = X.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    //if (is_gpu_backend() || use_gpu) {
    //    dispatch::cuda::ctzRcpp(result, X);
    //} else {
        dispatch::cpu::ctzRcpp(result, X);
    //}
    return result;
}