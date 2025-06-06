#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector Omega(const Rcpp::NumericMatrix & W);
    Rcpp::NumericVector omega(Rcpp::NumericMatrix & W);
}