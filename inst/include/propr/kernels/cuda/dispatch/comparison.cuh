#pragma once

namespace propr {
    namespace dispatch {
        namespace cuda {
            int count_less_than         (Rcpp::NumericVector x, double cutoff,propr::propr_context context);
            int count_greater_than      (Rcpp::NumericVector x, double cutoff,propr::propr_context context);
            int count_less_equal_than   (Rcpp::NumericVector x, double cutoff,propr::propr_context context);
            int count_greater_equal_than(Rcpp::NumericVector x, double cutoff,propr::propr_context context);
        }
    }
}