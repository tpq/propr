#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            int count_less_than(Rcpp::NumericVector& x, double cutoff, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            int count_greater_than(Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            int count_less_equal_than(Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            int count_greater_equal_than(Rcpp::NumericVector& x, double cutoff,propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);

            void* count_values_beyond_thresholds_begin(Rcpp::NumericVector& cutoffs, bool direct, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void count_values_beyond_thresholds_accumulate(void* counter, Rcpp::NumericVector& values, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            Rcpp::NumericVector count_values_beyond_thresholds_end(void* counter, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void count_values_beyond_thresholds_destroy(void* counter);
        }
    }
}
