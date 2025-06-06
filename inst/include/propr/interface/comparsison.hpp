#pragma once

#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    int count_less_than         (Rcpp::NumericVector x, double cutoff);
    int count_greater_than      (Rcpp::NumericVector x, double cutoff);
    int count_less_equal_than   (Rcpp::NumericVector x, double cutoff);
    int count_greater_equal_than(Rcpp::NumericVector x, double cutoff);
}
