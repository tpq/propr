#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector getOR    (const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, bool use_gpu = false);
    Rcpp::NumericVector getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm, bool use_gpu = false);
    Rcpp::NumericMatrix permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p = 100, bool use_gpu = false);
    Rcpp::List getFDR            (double actual, const Rcpp::NumericVector& permuted, bool use_gpu = false);
    Rcpp::IntegerMatrix getG     (const Rcpp::IntegerVector& Gk, bool use_gpu = false);
    Rcpp::NumericVector graflex  (const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p = 100, bool use_gpu = false);
}