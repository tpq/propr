#include <Rcpp.h>
using namespace Rcpp;

// Count the number of elements in vector `x` that are strictly less than the
// `cutoff`.
//
// [[Rcpp::export]]
int count_less_than(NumericVector x, double cutoff) {
  // See `?Memory-limits`: R vectors are limited in size to the max value of a
  // signed int, as they use signed int for their size.  If all values of `x`
  // are greater than cutoff, then we would get a count of INT_MAX, so the count
  // shouldn't overflow.
  int count = 0;
  int len = x.size();

  for (int i = 0; i < len; ++i) {
    // Returns 1 if it's less than cutoff, zero otherwise.  Add it to the count.
    count += x[i] < cutoff;
  }

  return count;
}

// Count the number of elements in vector `x` that are strictly greater than the
// `cutoff`.
//
// [[Rcpp::export]]
int count_greater_than(NumericVector x, double cutoff) {
  int count = 0;
  int len = x.size();

  for (int i = 0; i < len; ++i) {
    // Returns 1 if it's less than cutoff, zero otherwise.  Add it to the count.
    count += x[i] > cutoff;
  }

  return count;
}
