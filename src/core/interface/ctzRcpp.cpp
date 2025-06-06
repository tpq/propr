// src/core/interface/ctzRcpp.cpp
#include <Rcpp.h>
#include "../../include/interface/ctzRcpp.hpp"
#include "../../include/interface/device_selector.hpp"
#include "../../include/context.h"

// Include CPU dispatch headers
#include "../../include/kernels/cpu/dispatch/ctzRcpp.hpp"

// Include CUDA dispatch headers
#include "../../include/kernels/cuda/dispatch/ctzRcpp.cuh"


// [[Rcpp::export]]
Rcpp::NumericVector propr::ctzRcpp(Rcpp::NumericMatrix & X) {
    int nfeats = X.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result(llt);

    if (propr::is_gpu_backend()) {
        propr::propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed for ctzRcpp: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::ctzRcpp(X, result);
        } else {
            context.stream = stream;
            propr::dispatch::cuda::ctzRcpp(X, result, context);
            cudaStreamDestroy(stream);
        }
    } else {
        propr::dispatch::cpu::ctzRcpp(X, result);
    }
    return result;
}