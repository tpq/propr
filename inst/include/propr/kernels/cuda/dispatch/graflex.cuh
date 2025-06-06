#pragma once

#include <Rcpp.h>
#include <propr/context.h>


namespace propr {
    namespace dispatch {
        namespace cuda {
            void getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, Rcpp::NumericVector& out, propr::propr_context context);
            void getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm, Rcpp::NumericVector&out, propr::propr_context context);
            void permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p, Rcpp::NumericMatrix& out, propr::propr_context context);
            void getFDR(double actual, const Rcpp::NumericVector& permuted, Rcpp::List& out, propr::propr_context context);
            void getG(const Rcpp::IntegerVector& Gk, Rcpp::IntegerMatrix& out, propr::propr_context context);
            void graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p, Rcpp::NumericVector& out, propr::propr_context context) ;
        }
    }
}
