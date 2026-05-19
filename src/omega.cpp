#include <Rcpp.h>

#include <propr/interface/omega.hpp>
#include <propr/kernels/cpu/dispatch/omega.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector omega(Rcpp::NumericMatrix& W, Rcpp::String backend = "auto") {
    const std::size_t nfeats = W.ncol();
    const std::size_t llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::dof_global(result, W);
    } else {
        dispatch::cpu::dof_global(result, W);
    }

    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector Omega(Rcpp::NumericMatrix& W, Rcpp::String backend = "auto") {
    const std::size_t nfeats = W.ncol();
    const std::size_t llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::dof_population(result, W);
    } else {
        dispatch::cpu::dof_population(result, W);
    }

    return result;
}
