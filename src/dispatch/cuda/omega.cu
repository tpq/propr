#include <Rcpp.h>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/cuda_checks.h>

#include <propr/context.h>
#include <propr/kernels/cuda/dispatch/omega.cuh>

using namespace Rcpp;
using namespace propr;

void dispatch::cuda::dof_global(NumericVector& out, const NumericMatrix & W, propr_context context) {
    int nfeats = W.ncol();
    size_t llt = size_t(nfeats) * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("dof_global is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::dof_population(NumericVector& out, const NumericMatrix & W, propr_context context) {
    int nfeats = W.ncol();
    size_t llt = size_t(nfeats) * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("dof_population is not implemented in CUDA. Falling back to CPU via dispatcher.");
}