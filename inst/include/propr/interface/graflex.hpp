#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericVector getOR    (const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G);
    Rcpp::NumericVector getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm);
    Rcpp::NumericMatrix permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p = 100);
    Rcpp::List getFDR            (double actual, const Rcpp::NumericVector& permuted);
    Rcpp::IntegerMatrix getG     (const Rcpp::IntegerVector& Gk);
    Rcpp::NumericVector graflex  (const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p = 100);
}