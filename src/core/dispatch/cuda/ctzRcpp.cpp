#include <Rcpp.h>
#include "../../../include/kernels/cuda/dispatch/ctzRcpp.cuh"

#include "utils.hpp"

using namespace Rcpp;
using namespace propr;

void dispatch::cuda::ctzRcpp(NumericMatrix & X, NumericVector& out, propr::propr_context context){
    Rcpp::stop("Not implemented in CUDA for ctzRcpp.");
    int nfeats = X.ncol();
    int nsubjs = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    Rcpp::NumericVector zeroes(nfeats); // Count zero frequency per feature
    Rcpp::NumericVector result(llt);    // Count joint zero frequency
    return result;
}