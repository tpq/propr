#include <Rcpp.h>
#include "../../../include/kernels/cpu/dispatch/comparison.hpp"

using namespace propr;

int dispatch::cpu::count_less_than(Rcpp::NumericVector &x, double cutoff) {
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] < cutoff;
  }
  return count;
}

int dispatch::cpu::count_greater_than(Rcpp::NumericVector &x, double cutoff) {
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] > cutoff;
  }
  return count;
}

int dispatch::cpu::count_less_equal_than(Rcpp::NumericVector &x, double cutoff) {
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] <= cutoff;
  }
  return count;
}

int dispatch::cpu::count_greater_equal_than(Rcpp::NumericVector &x, double cutoff) {
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] >= cutoff;
  }
  return count;
}