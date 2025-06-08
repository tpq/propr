#include <Rcpp.h>
#include <propr/interface/graflex.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/graflex.hpp>
#include <propr/kernels/cuda/dispatch/graflex.cuh>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G) {
    Rcpp::NumericVector result(8); // Assuming getOR returns a vector of size 8
     if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getOR: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::getOR(A, G, result);
        } else {
            context.stream = stream;
            dispatch::cuda::getOR(A, G,result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::getOR(A, G, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm) {
    Rcpp::NumericVector result(8);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getORperm: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::getORperm(A, G, perm, result);
        } else {
            context.stream = stream;
            dispatch::cuda::getORperm(A, G, perm, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::getORperm(A, G, perm, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p) {
     if (is_gpu_backend()) {
        propr_context context;
        Rcpp::NumericMatrix result(p, 8);
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for permuteOR: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::permuteOR(A, G, p, result);
        } else {
            context.stream = stream;
             dispatch::cuda::permuteOR(A, G, p, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        Rcpp::NumericMatrix result(p, 8);
        dispatch::cpu::permuteOR(A, G, p, result);
        return result;
    }
}

// [[Rcpp::export]]
Rcpp::List getFDR(double actual, const Rcpp::NumericVector& permuted) {
    Rcpp::List result;
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getFDR: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::getFDR(actual, permuted, result);
        } else {
            context.stream = stream;
            dispatch::cuda::getFDR(actual, permuted, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::getFDR(actual, permuted, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getG(const Rcpp::IntegerVector& Gk) {
    int n = Gk.size();
    Rcpp::IntegerMatrix result(n, n);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for getG: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::getG(Gk, result);
        } else {
            context.stream = stream;
            dispatch::cuda::getG(Gk, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::getG(Gk, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p) {
    Rcpp::NumericVector result(8);
    if (is_gpu_backend()) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for graflex: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::graflex(A, Gk, p, result);
        } else {
            context.stream = stream;
            dispatch::cuda::graflex(A, Gk, p, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        dispatch::cpu::graflex(A, Gk, p, result);
    }
    return result;
}