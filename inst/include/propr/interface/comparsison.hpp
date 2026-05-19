#pragma once

#include <Rcpp.h>

namespace propr {
    int count_less_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_greater_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_less_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_greater_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
}
