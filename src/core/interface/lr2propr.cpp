#include <Rcpp.h>
#include "../../include/interface/lr2propr.hpp"
#include "../../include/interface/device_selector.hpp"
#include "../../include/context.h"

#include "../../include/kernels/cpu/dispatch/lr2propr.hpp"
#include "../../include/kernels/cuda/dispatch/lr2propr.cuh"


// [[Rcpp::export]]
Rcpp::NumericMatrix propr::lr2vlr(Rcpp::NumericMatrix lr) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for lr2vlr: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::lr2vlr(lr, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::lr2vlr(lr, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::lr2vlr(lr, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix propr::lr2phi(Rcpp::NumericMatrix lr) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for lr2phi: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::lr2phi(lr, result);
        } else {
            context.stream = stream;
            result = propr::dispatch::cuda::lr2phi(lr,result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::lr2phi(lr, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix propr::lr2rho(Rcpp::NumericMatrix lr) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for lr2rho: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::lr2rho(lr, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::lr2rho(lr, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::lr2rho(lr, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix propr::lr2phs(Rcpp::NumericMatrix lr) {
    int nfeats = lr.ncol();
    Rcpp::NumericMatrix result(nfeats, nfeats);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for lr2phs: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::lr2phs(lr, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::lr2phs(lr, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::lr2phs(lr, result);
    }
    return result;
}