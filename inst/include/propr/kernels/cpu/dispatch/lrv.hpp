#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu   {
            void lrv(Rcpp::NumericVector& out,
                    Rcpp::NumericMatrix &Y,
                    Rcpp::NumericMatrix &W,
                    bool weighted = false,
                    double a = NA_REAL,
                    Rcpp::NumericMatrix Yfull = Rcpp::NumericMatrix(1, 1).nrow(),
                    Rcpp::NumericMatrix Wfull = Rcpp::NumericMatrix(1, 1).nrow());
        }
    }
}