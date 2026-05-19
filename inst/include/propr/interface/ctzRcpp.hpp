#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector ctzRcpp(Rcpp::NumericMatrix& X, Rcpp::String backend = "auto");
}
