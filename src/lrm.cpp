#include <Rcpp.h>

#include <propr/interface/lrm.hpp>
#include <propr/kernels/cpu/dispatch/lrm.hpp>
#include <propr/runtime/cuda_executor.hpp>
#include <propr/runtime/dispatch.hpp>

using namespace propr;

// [[Rcpp::export]]
Rcpp::NumericVector lrm(
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
                runtime::cuda_executor::lrm_alpha_weighted(result_vec, Y, W, a, Yfull, Wfull);
            } else {
                runtime::cuda_executor::lrm_alpha(result_vec, Y, a, Yfull);
            }
        } else {
            if (weighted) {
                runtime::cuda_executor::lrm_weighted(result_vec, Y, Wfull);
            } else {
                runtime::cuda_executor::lrm_basic(result_vec, Y);
            }
        }
    } else {
        dispatch::cpu::lrm(result_vec, Y, W, weighted, a, Yfull, Wfull);
    }

    return result_vec;
}
