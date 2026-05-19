#include <Rcpp.h>

#include <propr/interface/ctzRcpp.hpp>
#include <propr/kernels/cpu/dispatch/ctzRcpp.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector ctzRcpp(Rcpp::NumericMatrix& X, Rcpp::String backend = "auto") {
    const int nfeats = X.ncol();
    const int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::ctzRcpp(result, X);
    } else {
        dispatch::cpu::ctzRcpp(result, X);
    }
    return result;
}
