#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void getOR(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G);
            void getORperm(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm);
            void permuteOR(Rcpp::NumericMatrix& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p);
            void getFDR(Rcpp::List& out, double actual, const Rcpp::NumericVector& permuted);
            void getG(Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& Gk);
            void graflex(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p);
        }
    }
}