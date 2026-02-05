#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector ctzRcpp(Rcpp::NumericMatrix & X, bool use_gpu = false);
}