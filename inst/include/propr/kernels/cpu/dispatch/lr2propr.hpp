#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void lr2vlr(Rcpp::NumericMatrix &lr, Rcpp::NumericMatrix& out);
            void lr2phi(Rcpp::NumericMatrix &lr, Rcpp::NumericMatrix& out);
            void lr2rho(Rcpp::NumericMatrix &lr, Rcpp::NumericMatrix& out);
            void lr2phs(Rcpp::NumericMatrix &lr, Rcpp::NumericMatrix& out);
        }
    }
}