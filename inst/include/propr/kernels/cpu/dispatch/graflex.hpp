#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            
            __forceinline__
            void getOR(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G);
            
            __forceinline__
            void getORperm(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm);
            
            __forceinline__
            void permuteOR(Rcpp::NumericMatrix& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p);
            
            __forceinline__
            void getFDR(Rcpp::List& out, double actual, const Rcpp::NumericVector& permuted);
            
            __forceinline__
            void getG(Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& Gk);
            
            __forceinline__
            void graflex(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p);
        }
    }
}