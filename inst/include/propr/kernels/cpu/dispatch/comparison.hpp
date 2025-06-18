#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
      namespace cpu {
          
        __forceinline__
        int count_less_than         (Rcpp::NumericVector &x, double cutoff);
          
        __forceinline__  
        int count_greater_than      (Rcpp::NumericVector &x, double cutoff);
          
        __forceinline__
        int count_less_equal_than   (Rcpp::NumericVector &x, double cutoff);
          
        __forceinline__
        int count_greater_equal_than(Rcpp::NumericVector &x, double cutoff);
      }
    }
}