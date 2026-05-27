#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <propr/kernels/cpu/dispatch/comparison.hpp>
#include <propr/utils/profilers/host_exclusive_profiler.hpp>

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

// begin: these types are kind of shared between the cpu and the gpu, we need to restructure things more appropriately
using threshold_count_t = unsigned long long;

struct threshold_cutoff {
  double value;
  R_xlen_t index;
};

struct threshold_group {
  std::vector<threshold_cutoff> cutoffs;
  std::vector<double> values;
  std::vector<threshold_count_t> buckets;
};

struct cpu_threshold_counter {
  R_xlen_t ncutoffs;
  threshold_group less;
  threshold_group greater;
};
// end:

static void prepare_group(threshold_group& group) {
  std::stable_sort( group.cutoffs.begin(), group.cutoffs.end(),
    [](const threshold_cutoff& a, const threshold_cutoff& b) {
      return a.value < b.value;
    });
  group.values.resize(group.cutoffs.size());
  group.buckets.assign(group.cutoffs.size() + 1, 0);
  for (size_t i = 0; i < group.cutoffs.size(); ++i) {
    group.values[i] = group.cutoffs[i].value;
  }
}

void* dispatch::cpu::count_values_beyond_thresholds_begin( Rcpp::NumericVector& cutoffs, bool direct) {
  auto* counter = new cpu_threshold_counter{};
  counter->ncutoffs = cutoffs.size();

  for (R_xlen_t i = 0; i < cutoffs.size(); ++i) {
    const double cutoff = cutoffs[i];
    if (std::isnan(cutoff))  continue;
    if (direct && cutoff >= 0.0) {
      counter->greater.cutoffs.push_back({cutoff, i});
    } else {
      counter->less.cutoffs.push_back({cutoff, i});
    }
  }

  prepare_group(counter->less);
  prepare_group(counter->greater);
  return counter;
}

void dispatch::cpu::count_values_beyond_thresholds_accumulate( void* raw_counter, Rcpp::NumericVector& values) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel");

  auto* counter = static_cast<cpu_threshold_counter*>(raw_counter);

  if (!counter->less.cutoffs.empty()) {
    for (R_xlen_t i = 0; i < values.size(); ++i) {
      const double value = values[i];
      if (std::isnan(value)) continue;
      const auto it = std::upper_bound( counter->less.values.begin(), counter->less.values.end(), value);
      ++counter->less.buckets[static_cast<size_t>(it - counter->less.values.begin())];
    }
  }

  if (!counter->greater.cutoffs.empty()) {
    for (R_xlen_t i = 0; i < values.size(); ++i) {
      const double value = values[i];
      if (std::isnan(value))  continue;
      const auto it = std::lower_bound( counter->greater.values.begin(), counter->greater.values.end(), value);
      ++counter->greater.buckets[ static_cast<size_t>(it - counter->greater.values.begin())];
    }
  }
}

Rcpp::NumericVector dispatch::cpu::count_values_beyond_thresholds_end( void* raw_counter) {
  auto* counter = static_cast<cpu_threshold_counter*>(raw_counter);
  Rcpp::NumericVector out(counter->ncutoffs, NA_REAL);

  threshold_count_t prefix = 0;
  for (size_t i = 0; i < counter->less.cutoffs.size(); ++i) {
    prefix += counter->less.buckets[i];
    out[counter->less.cutoffs[i].index] = static_cast<double>(prefix);
  }

  threshold_count_t suffix = 0;
  for (size_t i = counter->greater.cutoffs.size(); i > 0; --i) {
    suffix += counter->greater.buckets[i];
    out[counter->greater.cutoffs[i - 1].index] = static_cast<double>(suffix);
  }

  return out;
}

void dispatch::cpu::count_values_beyond_thresholds_destroy(void* raw_counter) {
  delete static_cast<cpu_threshold_counter*>(raw_counter);
}
