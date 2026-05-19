#include <Rcpp.h>

#include <propr/interface/graflex.hpp>
#include <propr/kernels/cpu/dispatch/graflex.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, Rcpp::String backend = "auto") {
    Rcpp::NumericVector result(8);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::getOR(result, A, G);
    } else {
        dispatch::cpu::getOR(result, A, G);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector getORperm(
    const Rcpp::IntegerMatrix& A,
    const Rcpp::IntegerMatrix& G,
    const Rcpp::IntegerVector& perm,
    Rcpp::String backend = "auto") {
    Rcpp::NumericVector result(8);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::getORperm(result, A, G, perm);
    } else {
        dispatch::cpu::getORperm(result, A, G, perm);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p, Rcpp::String backend = "auto") {
    Rcpp::NumericMatrix result(p, 8);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::permuteOR(result, A, G, p);
    } else {
        dispatch::cpu::permuteOR(result, A, G, p);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List getFDR(double actual, const Rcpp::NumericVector& permuted, Rcpp::String backend = "auto") {
    Rcpp::List result;
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::getFDR(result, actual, permuted);
    } else {
        dispatch::cpu::getFDR(result, actual, permuted);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getG(const Rcpp::IntegerVector& Gk, Rcpp::String backend = "auto") {
    const int n = Gk.size();
    Rcpp::IntegerMatrix result(n, n);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::getG(result, Gk);
    } else {
        dispatch::cpu::getG(result, Gk);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p, Rcpp::String backend = "auto") {
    Rcpp::NumericVector result(8);
    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        runtime::cuda_executor::graflex(result, A, Gk, p);
    } else {
        dispatch::cpu::graflex(result, A, Gk, p);
    }
    return result;
}
