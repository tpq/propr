#include <Rcpp.h>

namespace propr {
    namespace dispatch {
      namespace cpu {
          int count_less_than         (Rcpp::NumericVector &x, double cutoff);
          int count_greater_than      (Rcpp::NumericVector &x, double cutoff);
          int count_less_equal_than   (Rcpp::NumericVector &x, double cutoff);
          int count_greater_equal_than(Rcpp::NumericVector &x, double cutoff);
      }
    }
}