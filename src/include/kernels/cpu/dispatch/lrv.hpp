#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu   {
            void lrv(Rcpp::NumericMatrix &Y,
                  Rcpp::NumericMatrix &W,
                  bool weighted = false,
                  double a = NA_REAL,
                  Rcpp::NumericMatrix& Yfull,
                  Rcpp::NumericMatrix& Wfull,
                  Rcpp::NumericVector& out);
        }
    }
}