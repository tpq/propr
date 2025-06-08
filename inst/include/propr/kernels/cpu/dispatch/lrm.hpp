#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void lrm(Rcpp::NumericVector& out,
                    Rcpp::NumericMatrix &Y, 
                    Rcpp::NumericMatrix &W,
                    bool weighted, double a,
                    Rcpp::NumericMatrix Yfull =  Rcpp::NumericMatrix(1, 1).nrow(),
                    Rcpp::NumericMatrix Wfull =  Rcpp::NumericMatrix(1, 1).nrow());
        }
    }
}