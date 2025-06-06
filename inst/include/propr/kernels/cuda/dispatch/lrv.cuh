#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            void lrv_basic(
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericVector &out,
                propr_context context);

            void lrv_weighted(
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                Rcpp::NumericVector& out,
                propr_context context);

            void lrv_alpha(
                Rcpp::NumericMatrix &Y,
                const double a,
                Rcpp::NumericMatrix& Yfull,
                Rcpp::NumericVector &out,
                propr_context context);

            void lrv_alpha_weighted(
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                const double a,
                Rcpp::NumericMatrix& Yfull,
                Rcpp::NumericMatrix& Wfull,
                Rcpp::NumericVector& out, propr_context context);

        }
    } 
}
