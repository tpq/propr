#include <Rcpp.h>
#include "../../../include/kernels/cuda/dispatch/comparison.cuh"
#include "context.h"

using namespace Rcpp;
using namespace propr;

int dispatch::cuda::count_less_than(NumericVector x, double cutoff, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for count_less_than.");
    return 0;
}

int dispatch::cuda::count_greater_than(NumericVector x, double cutoff, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for count_greater_than.");
    return 0;
}

int dispatch::cuda::count_less_equal_than(NumericVector x, double cutoff, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for count_less_equal_than.");
    return 0;
}

int dispatch::cuda::count_greater_equal_than(NumericVector x, double cutoff, propr::propr_context context) {
    Rcpp::stop("Not implemented in CUDA for count_greater_equal_than.");
    return 0;
}