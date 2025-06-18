#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {

            __forceinline__
            void lrv_basic(
                Rcpp::NumericVector &out,
                Rcpp::NumericMatrix &Y,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lrv_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lrv_alpha(
                Rcpp::NumericVector &out,
                Rcpp::NumericMatrix &Y,
                const double a,
                Rcpp::NumericMatrix& Yfull,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lrv_alpha_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                const double a,
                Rcpp::NumericMatrix& Yfull,
                Rcpp::NumericMatrix& Wfull,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);

        }
    } 
}
