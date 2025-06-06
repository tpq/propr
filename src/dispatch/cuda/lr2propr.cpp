#include <Rcpp.h>
#include <math.h>
#include <propr/kernels/cuda/dispatch/lr2propr.cuh>

#include <propr/context.h>

using namespace Rcpp;
using namespace propr;

void dispatch::cuda::lr2vlr(NumericMatrix& lr, Rcpp::NumericMatrix &out, propr_context context) {
    Rcpp::stop("Not implemented in CUDA for lr2vlr.");
    // Calculate variation matrix
    NumericMatrix x = clone(lr);
}

void dispatch::cuda::lr2phi(NumericMatrix& lr, Rcpp::NumericMatrix &out, propr_context context) {
    Rcpp::stop("Not implemented in CUDA for lr2phi.");
    // Make vlr from log-ratio data
    NumericMatrix x = clone(lr);
    // NumericMatrix mat = lr2vlr(x);
}

void dispatch::cuda::lr2rho(NumericMatrix& lr,  Rcpp::NumericMatrix &out,propr_context context) {
    Rcpp::stop("Not implemented in CUDA for lr2rho.");
    NumericMatrix x = clone(lr);
    // NumericMatrix mat = lr2vlr(x);
}

void dispatch::cuda::lr2phs(NumericMatrix &lr, Rcpp::NumericMatrix &out, propr_context context) {
    Rcpp::stop("Not implemented in CUDA for lr2phs.");
    // Calculate phs = (1 - rho) / (1 + rho)
    // NumericMatrix mat = lr2rho(lr); // This call will now also stop if lr2rho is not implemented.
}