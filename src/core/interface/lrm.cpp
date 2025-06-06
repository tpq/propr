#include <Rcpp.h>

#include "../../include/interface/lrm.hpp"
#include "../../include/interface/device_selector.hpp"
#include "../../include/context.h"

#include "../../include/kernels/cpu/dispatch/lrm.hpp"
#include "../../include/kernels/cuda/dispatch/lrm.cuh"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lrm(NumericMatrix &Y,
                  NumericMatrix &W,
                  bool weighted,
                  double a,
                  NumericMatrix Yfull,
                  NumericMatrix Wfull) {

    const bool use_gpu = propr::is_gpu_backend();

    int nfeats = Y.ncol();
    int N_pairs = nfeats * (nfeats - 1) / 2;
    NumericVector result_vec(N_pairs);

    if (use_gpu) {
        propr::propr_context context;
        cudaStream_t stream;
        const cudaError_t err = cudaStreamCreate(&stream);
        if (err != cudaSuccess) {
            Rcpp::warning("CUDA stream creation failed: %s. Falling back to CPU.", cudaGetErrorString(err));
            propr::dispatch::cpu::lrm(Y, W, weighted, a, Yfull, Wfull, result_vec);
            return result_vec;
        }
        context.stream = stream;

        if (!R_IsNA(a)) {
            if (weighted) {
                 propr::dispatch::cuda::lrm_alpha_weighted(Y, W, a, Yfull, Wfull, result_vec, context);
            } else {
                propr::dispatch::cuda::lrm_alpha(Y, a, Yfull, result_vec, context);
            }
        } else {
            if (weighted) {
                propr::dispatch::cuda::lrm_weighted(Y, Wfull, result_vec, context);
            } else {
                propr::dispatch::cuda::lrm_basic(Y, result_vec, context);
            }
        }
        cudaStreamDestroy(stream);
        return result_vec;

    } else {
        propr::dispatch::cpu::lrm(Y, W, weighted, a, Yfull, Wfull, result_vec);
        return result_vec;
    }
}