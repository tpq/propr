#include <Rcpp.h>
#include <propr/interface/ctzRcpp.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/ctzRcpp.hpp>
#include <propr/kernels/cuda/dispatch/ctzRcpp.cuh>


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