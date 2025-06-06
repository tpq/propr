#include <Rcpp.h>
#include "../../../include/kernels/cpu/dispatch/graflex.hpp"
#include "../../../include/utils.hpp"

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::getOR(const IntegerMatrix& A, const IntegerMatrix& G, NumericVector& out) {
  CHECK_VECTOR_SIZE(out, 8);
  int ncol = A.ncol();
  int a = 0, b = 0, c = 0, d = 0;

  for (int j = 0; j < ncol - 1; ++j) {
    for (int i = j + 1; i < ncol; ++i) {
      int a_val = A(i, j), g_val = G(i, j);
      a += (1 - a_val) * (1 - g_val);
      b += (1 - a_val) * g_val;
      c += a_val * (1 - g_val);
      d += a_val * g_val;
    }
  }

  double odds_ratio = static_cast<double>(a * d) / (b * c);
  double log_odds_ratio = std::log(odds_ratio);

  out[0] = a;
  out[1] = b;
  out[2] = c;
  out[3] = d;
  out[4] = odds_ratio;
  out[5] = log_odds_ratio;
  out[6] = R_NaN;
  out[7] = R_NaN;
}

void
dispatch::cpu::getORperm(const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm, NumericVector& out) {
  CHECK_VECTOR_SIZE(out, 8); // Check if output vector has correct size

  int ncol = A.ncol();
  int a = 0, b = 0, c = 0, d = 0;

  for (int j = 0; j < ncol - 1; ++j) {
    int pj = perm[j];
    for (int i = j + 1; i < ncol; ++i) {
      int a_val = A(perm[i], pj), g_val = G(i, j);
      a += (1 - a_val) * (1 - g_val);
      b += (1 - a_val) * g_val;
      c += a_val * (1 - g_val);
      d += a_val * g_val;
    }
  }

  double odds_ratio = static_cast<double>(a * d) / (b * c);
  double log_odds_ratio = std::log(odds_ratio);

  out[0] = a;
  out[1] = b;
  out[2] = c;
  out[3] = d;
  out[4] = odds_ratio;
  out[5] = log_odds_ratio;
  out[6] = R_NaN;
  out[7] = R_NaN;
}

void
dispatch::cpu::permuteOR(const IntegerMatrix& A, const IntegerMatrix& G, int p, NumericMatrix& out) {
  CHECK_MATRIX_DIMS(out, p, 8);
  int ncol = A.ncol();
  NumericVector or_tmp(8);
  for (int i = 0; i < p; ++i) {
    IntegerVector perm = sample(ncol, ncol, false) - 1;
    dispatch::cpu::getORperm(A, G, perm, or_tmp);
    for (int j = 0; j < 8; ++j) {
        out(i, j) = or_tmp[j];
    }
  }
}

void
dispatch::cpu::getFDR(double actual, const NumericVector& permuted, List& out) {
  int n = permuted.size();
  int count_over = 0;
  int count_under = 0;

  for (int i = 0; i < n; ++i) {
    double current = permuted[i];
    if (current >= actual) ++count_over;
    if (current <= actual) ++count_under;
  }

  double fdr_over = static_cast<double>(count_over) / n;
  double fdr_under = static_cast<double>(count_under) / n;

  out["over"] = fdr_over;
  out["under"] = fdr_under;
}

void
dispatch::cpu::getG(const IntegerVector& Gk, IntegerMatrix& out) {
  int n = Gk.size();
  CHECK_MATRIX_DIMS(out, n, n);
  for (int i = 0; i < n; ++i) {
    int gi = Gk[i];
    out(i, i) = gi * gi;
    for (int j = 0; j < i; ++j) {
      int value = gi * Gk[j];
      out(i, j) = value;
      out(j, i) = value;
    }
  }
}

void
dispatch::cpu::graflex(const IntegerMatrix& A, const IntegerVector& Gk, int p, NumericVector& out) {
  IntegerMatrix G_tmp(Gk.size(), Gk.size());
  dispatch::cpu::getG(Gk, G_tmp);

  NumericVector actual_tmp(8);
  dispatch::cpu::getOR(A, G_tmp, actual_tmp);

  CHECK_VECTOR_SIZE(out, actual_tmp.length());
  for (int i = 0; i < actual_tmp.length(); ++i) {
      out[i] = actual_tmp[i];
  }

  if (!std::isnan(actual_tmp(4))) {
    NumericMatrix permuted_tmp(p, 8);
    dispatch::cpu::permuteOR(A, G_tmp, p, permuted_tmp);

    List fdr_tmp; // Temporary list for FDR results
    dispatch::cpu::getFDR(actual_tmp(4), permuted_tmp(_, 4), fdr_tmp);
    out(6) = as<double>(fdr_tmp["under"]);
    out(7) = as<double>(fdr_tmp["over"]);
  }
}