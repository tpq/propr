#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {

            __forceinline__
            void lr2vlr(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &lr);
            
            __forceinline__
            void lr2phi(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &lr);
            
            __forceinline__
            void lr2rho(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &lr);
            
            __forceinline__
            void lr2phs(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &lr);
        }
    }
}