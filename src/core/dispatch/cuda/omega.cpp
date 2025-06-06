#include <Rcpp.h>
#include "../../../include/utils.hpp"
#include "../../../include/context.h"
#include "../../../include/kernels/cuda/dispatch/omega.cuh"

using namespace Rcpp;
using namespace propr;

void dispatch::cuda::dof_global(const NumericMatrix & W, NumericVector& out, propr_context context) {
    int nfeats = W.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("dof_global is not implemented in CUDA. Falling back to CPU via dispatcher.");
}

void dispatch::cuda::dof_population(const NumericMatrix & W, NumericVector& out, propr_context context) {
    int nfeats = W.ncol();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::stop("dof_population is not implemented in CUDA. Falling back to CPU via dispatcher.");
}