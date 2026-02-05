#include <Rcpp.h>
#include <propr/interface/graflex.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/graflex.hpp>
#include <propr/kernels/cuda/dispatch/graflex.cuh>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, bool use_gpu) {
    Rcpp::NumericVector result(8);
     if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::getOR(result, A, G);
    } else {
        dispatch::cpu::getOR(result, A, G);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm, bool use_gpu) {
    Rcpp::NumericVector result(8);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::getORperm(result, A, G, perm);
    } else {
        dispatch::cpu::getORperm(result, A, G, perm);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p, bool use_gpu) {
    Rcpp::NumericMatrix result(p, 8);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::permuteOR(result, A, G, p);
    } else {
        dispatch::cpu::permuteOR(result, A, G, p);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List getFDR(double actual, const Rcpp::NumericVector& permuted, bool use_gpu) {
    Rcpp::List result;
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::getFDR(result, actual, permuted);
    } else {
        dispatch::cpu::getFDR(result, actual, permuted);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getG(const Rcpp::IntegerVector& Gk, bool use_gpu) {
    int n = Gk.size();
    Rcpp::IntegerMatrix result(n, n);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::getG(result, Gk);
    } else {
        dispatch::cpu::getG(result, Gk);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p, bool use_gpu) {
    Rcpp::NumericVector result(8);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::graflex(result, A, Gk, p);
    } else {
        dispatch::cpu::graflex(result, A, Gk, p);
    }
    return result;
}