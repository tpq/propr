#pragma once

#include <Rcpp.h>

namespace propr {
    int count_less_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_greater_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_less_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");
    int count_greater_equal_than(Rcpp::NumericVector x, double cutoff, Rcpp::String backend = "auto");

    SEXP count_values_beyond_thresholds_begin(Rcpp::NumericVector cutoffs, bool direct, Rcpp::String backend = "auto");
    void count_values_beyond_thresholds_accumulate(SEXP counter, Rcpp::NumericVector values);
    Rcpp::NumericVector count_values_beyond_thresholds_end(SEXP counter);
}
