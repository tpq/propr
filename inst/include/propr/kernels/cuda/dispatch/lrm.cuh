#pragma once
#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            
            void lrm_basic(Rcpp::NumericVector& out,
                           Rcpp::NumericMatrix &Y,
                           propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);

            void lrm_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);

            void lrm_alpha(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix &Y,
                double a,
                Rcpp::NumericMatrix &Yfull,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);

            void lrm_alpha_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W, double a,
                Rcpp::NumericMatrix &Yfull,
                Rcpp::NumericMatrix &Wfull,
                propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}
