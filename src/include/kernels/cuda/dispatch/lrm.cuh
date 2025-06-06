#pragma once

#include "context.h"

namespace propr {
    namespace dispatch {
        namespace cuda {
            void lrm_basic(Rcpp::NumericMatrix &Y,
                           Rcpp::NumericVector& out,
                           propr::propr_context context);

            void lrm_weighted(
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W,
                Rcpp::NumericVector& out,
                propr_context context);

            void lrm_alpha(
                Rcpp::NumericMatrix &Y,
                double a,
                Rcpp::NumericMatrix &Yfull,
                Rcpp::NumericVector& out,
                propr_context context);

            void lrm_alpha_weighted(
                Rcpp::NumericMatrix &Y,
                Rcpp::NumericMatrix &W, double a,
                Rcpp::NumericMatrix &Yfull,
                Rcpp::NumericMatrix &Wfull,
                Rcpp::NumericVector& out,
                propr_context context);
        }
    }
}
