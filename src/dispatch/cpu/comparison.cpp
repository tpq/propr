#include <Rcpp.h>
#include <propr/kernels/cpu/dispatch/comparison.hpp>
#include <propr/utils/host_exclusive_profiler.hpp>

using namespace propr;

int dispatch::cpu::count_less_than(Rcpp::NumericVector &x, double cutoff) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] < cutoff;
  }
  return count;
}

int dispatch::cpu::count_greater_than(Rcpp::NumericVector &x, double cutoff) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] > cutoff;
  }
  return count;
}

int dispatch::cpu::count_less_equal_than(Rcpp::NumericVector &x, double cutoff) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] <= cutoff;
  }
  return count;
}

int dispatch::cpu::count_greater_equal_than(Rcpp::NumericVector &x, double cutoff) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int count = 0;
  int len = x.size();
  for (int i = 0; i < len; ++i) {
    count += x[i] >= cutoff;
  }
  return count;
}