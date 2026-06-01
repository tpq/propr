#pragma once

#include <Rcpp.h>

namespace propr {
    namespace dispatch {
      namespace cpu {
          int count_less_than         (Rcpp::NumericVector &x, double cutoff);
          int count_greater_than      (Rcpp::NumericVector &x, double cutoff);
          int count_less_equal_than   (Rcpp::NumericVector &x, double cutoff);
          int count_greater_equal_than(Rcpp::NumericVector &x, double cutoff);

          void* count_values_beyond_thresholds_begin(Rcpp::NumericVector& cutoffs, bool direct);
          void count_values_beyond_thresholds_accumulate(void* counter, Rcpp::NumericVector& values);
          Rcpp::NumericVector count_values_beyond_thresholds_end(void* counter);
          void count_values_beyond_thresholds_destroy(void* counter);
      }
    }
}
