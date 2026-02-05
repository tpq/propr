#include <Rcpp.h>

#include <propr/interface/lrv.hpp>
#include <propr/kernels/cpu/dispatch/lrv.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector lrv(
    Rcpp::NumericMatrix& Y,
    Rcpp::NumericMatrix& W,
    bool weighted,
    double a,
    Rcpp::NumericMatrix Yfull,
    Rcpp::NumericMatrix Wfull,
    Rcpp::String backend = "auto") {
    const int nfeats = Y.ncol();
    const int n_pairs = nfeats * (nfeats - 1) / 2;
    Rcpp::NumericVector result_vec(n_pairs);

    if (runtime::resolve_backend(backend) == runtime::Backend::CUDA) {
        if (!R_IsNA(a)) {
            if (weighted) {
                runtime::cuda_executor::lrv_alpha_weighted(result_vec, Y, W, a, Yfull, Wfull);
            } else {
                runtime::cuda_executor::lrv_alpha(result_vec, Y, a, Yfull);
            }
        } else {
            if (weighted) {
                runtime::cuda_executor::lrv_weighted(result_vec, Y, W);
            } else {
                runtime::cuda_executor::lrv_basic(result_vec, Y);
            }
        }
    } else {
        dispatch::cpu::lrv(result_vec, Y, W, weighted, a, Yfull, Wfull);
    }

    return result_vec;
}
