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
        if (!R_IsNA(a)) { // Alpha-transformed
            if (weighted) {
                dispatch::cuda::lrv_alpha_weighted(result_vec, Y, W, a, Yfull, Wfull);
            } else {
                dispatch::cuda::lrv_alpha(result_vec, Y, a, Yfull);
            }
        } else { // Non-transformed (log)
            if (weighted) {
                dispatch::cuda::lrv_weighted(result_vec, Y, W);
            } else {
                dispatch::cuda::lrv_basic(result_vec, Y);
            }
        }
    } else {
        dispatch::cpu::lrv(result_vec, Y, W, weighted, a, Yfull, Wfull);
    }
    return result_vec;
}