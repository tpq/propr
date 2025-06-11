#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector Omega(const Rcpp::NumericMatrix & W, bool use_gpu = false);
    Rcpp::NumericVector omega(Rcpp::NumericMatrix & W, bool use_gpu = false);
}