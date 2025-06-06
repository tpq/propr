#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void lrm(Rcpp::NumericMatrix &Y, Rcpp::NumericMatrix &W,
                     bool weighted, double a,
                     Rcpp::NumericMatrix& Yfull,
                     Rcpp::NumericMatrix& Wfull,
                     Rcpp::NumericVector& out);
        }
    }
}