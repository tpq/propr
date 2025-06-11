#include <Rcpp.h>
#include <propr/interface/lr2propr.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/lr2propr.hpp>
#include <propr/kernels/cuda/dispatch/lr2propr.cuh>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2vlr(Rcpp::NumericMatrix lr, bool use_gpu) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::lr2vlr(result, lr);
    } else {
        dispatch::cpu::lr2vlr(result, lr);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2phi(Rcpp::NumericMatrix lr, bool use_gpu) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::lr2phi(result, lr);
    } else {
        dispatch::cpu::lr2phi(result, lr);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2rho(Rcpp::NumericMatrix lr, bool use_gpu) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);
    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::lr2rho(result, lr);
    } else {
        dispatch::cpu::lr2rho(result, lr);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix lr2phs(Rcpp::NumericMatrix lr, bool use_gpu) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (is_gpu_backend() || use_gpu) {
        dispatch::cuda::lr2phs(result, lr);
    } else {
        dispatch::cpu::lr2phs(result, lr);
    }
    return result;
}