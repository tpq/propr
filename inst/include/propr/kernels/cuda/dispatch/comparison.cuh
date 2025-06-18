#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            
            __forceinline__
            int count_less_than         (Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            int count_greater_than      (Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            int count_less_equal_than   (Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            int count_greater_equal_than(Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}