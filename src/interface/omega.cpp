#include <Rcpp.h>
#include <propr/interface/omega.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/omega.hpp>
#include <propr/kernels/cuda/dispatch/omega.cuh>


// [[Rcpp::export]]
Rcpp::NumericVector propr::Omega(const Rcpp::NumericMatrix & W) {
    int nfeats = W.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for Omega: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::dof_population(W, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::dof_population(W, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::dof_population(W, result);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector propr::omega(Rcpp::NumericMatrix & W) {
    int nfeats = W.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for omega: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::dof_global(W, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::dof_global(W, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::dof_global(W, result);
    }
    return result;
}