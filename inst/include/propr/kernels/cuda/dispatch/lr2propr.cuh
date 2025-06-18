#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            
            __forceinline__
            void lr2vlr(Rcpp::NumericMatrix &out, Rcpp::NumericMatrix& lr, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lr2phi(Rcpp::NumericMatrix &out, Rcpp::NumericMatrix& lr, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lr2rho(Rcpp::NumericMatrix &out, Rcpp::NumericMatrix &lr, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lr2phs(Rcpp::NumericMatrix &out, Rcpp::NumericMatrix& lr, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}
