#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector Omega(Rcpp::NumericMatrix& W, Rcpp::String backend = "auto");
    Rcpp::NumericVector omega(Rcpp::NumericMatrix& W, Rcpp::String backend = "auto");
}
