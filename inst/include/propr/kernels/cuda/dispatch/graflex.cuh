#pragma once

#include <Rcpp.h>
#include <propr/context.h>


namespace propr {
    namespace dispatch {
        namespace cuda {
            void getOR(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void getORperm(Rcpp::NumericVector&out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void permuteOR(Rcpp::NumericMatrix& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void getFDR(Rcpp::List& out, double actual, const Rcpp::NumericVector& permuted, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void getG(Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& Gk, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void graflex(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT) ;
        }
    }
}
