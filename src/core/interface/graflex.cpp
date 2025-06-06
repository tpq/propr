// src/core/interface/graflex.cpp
#include <Rcpp.h>
#include "../../include/interface/graflex.hpp"
#include "../../include/interface/device_selector.hpp"
#include "../../include/context.h"

// Include CPU dispatch headers
#include "../../include/kernels/cpu/dispatch/graflex.hpp"

// Include CUDA dispatch headers
#include "../../include/kernels/cuda/dispatch/graflex.cuh"


// [[Rcpp::export]]
Rcpp::NumericVector propr::getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G) {
    Rcpp::NumericVector result(8); // Assuming getOR returns a vector of size 8
     if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getOR: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::getOR(A, G, result);
        } else {
            context.stream = stream;
            result = propr::dispatch::cuda::getOR(A, G, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::getOR(A, G, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector propr::getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm) {
    Rcpp::NumericVector result(8);
    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getORperm: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::getORperm(A, G, perm, result);
        } else {
            context.stream = stream;
            // Assuming cuda::getORperm returns NumericVector
            result = propr::dispatch::cuda::getORperm(A, G, perm, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::getORperm(A, G, perm, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix propr::permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p) {
     if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for permuteOR: %s. Falling back to CPU.", cudaGetErrorString(err));
            Rcpp::NumericMatrix result(p, 8);
            propr::dispatch::cpu::permuteOR(A, G, p, result);
            return result;
        } else {
            context.stream = stream;
            Rcpp::NumericMatrix result = propr::dispatch::cuda::permuteOR(A, G, p, context);
            cudaStreamDestroy(stream);
            return result;
        }
    } else {
        Rcpp::NumericMatrix result(p, 8);
        propr::dispatch::cpu::permuteOR(A, G, p, result);
        return result;
    }
}

// [[Rcpp::export]]
Rcpp::List propr::getFDR(double actual, const Rcpp::NumericVector& permuted) {
    Rcpp::List result;
    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getFDR: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::getFDR(actual, permuted, result);
        } else {
            context.stream = stream;
            result = propr::dispatch::cuda::getFDR(actual, permuted, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::getFDR(actual, permuted, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix propr::getG(const Rcpp::IntegerVector& Gk) {
    int n = Gk.size();
    Rcpp::IntegerMatrix result(n, n);
    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getG: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::getG(Gk, result);
        } else {
            context.stream = stream;
            result = propr::dispatch::cuda::getG(Gk, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::getG(Gk, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector propr::graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p) {
    Rcpp::NumericVector result(8);
    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for graflex: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::graflex(A, Gk, p, result);
        } else {
            context.stream = stream;
            result = propr::dispatch::cuda::graflex(A, Gk, p, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::graflex(A, Gk, p, result);
    }
    return result;
}