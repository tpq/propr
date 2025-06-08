#include <Rcpp.h>
#include <propr/interface/lrv.hpp>
#include <propr/interface/device_selector.hpp>
#include <propr/context.h>

#include <propr/kernels/cpu/dispatch/lrv.hpp>
#include <propr/kernels/cuda/dispatch/lrv.cuh>

using namespace Rcpp;
using namespace propr;

// [[Rcpp::export]]
NumericVector lrv(NumericMatrix &Y,
                  NumericMatrix &W,
                  bool weighted,
                  double a,
                  NumericMatrix Yfull,
                  NumericMatrix Wfull) {

    bool use_gpu = is_gpu_backend();

    int nfeats = Y.ncol();
    int N_pairs = nfeats * (nfeats - 1) / 2;
    NumericVector result_vec(N_pairs);

    if (use_gpu) {
        propr_context context;
        cudaStream_t stream;
        cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed: %s. Falling back to CPU.", cudaGetErrorString(err));
            dispatch::cpu::lrv(Y, W, weighted, a, Yfull, Wfull, result_vec);
            return result_vec;
        }
        context.stream = stream;

        if (!R_IsNA(a)) { // Alpha-transformed
            if (weighted) {
                dispatch::cuda::lrv_alpha_weighted(Y, W, a, Yfull, Wfull, result_vec, context);
            } else {
                dispatch::cuda::lrv_alpha(Y, a, Yfull, result_vec, context);
            }
        } else { // Non-transformed (log)
            if (weighted) {
                dispatch::cuda::lrv_weighted(Y, W, result_vec, context);
            } else {
                dispatch::cuda::lrv_basic(Y, result_vec, context);
            }
        }
        cudaStreamDestroy(stream);
        return result_vec;
    } else {
        dispatch::cpu::lrv(Y, W, weighted, a, Yfull, Wfull, result_vec);
        return result_vec;
    }
}