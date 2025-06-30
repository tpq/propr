#include <Rcpp.h>
#include <propr/kernels/cpu/dispatch/graflex.hpp>
#include <propr/utils.hpp>

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::getOR(NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G) {
  CHECK_VECTOR_SIZE(out, 8);
  int ncol = A.ncol();
  int a = 0, b = 0, c = 0, d = 0;

  for (int j = 0; j < ncol - 1; ++j) {
    for (int i = j + 1; i < ncol; ++i) {
      int a_val = A(i, j);
      int g_val = G(i, j);
      a += (1 - a_val) * (1 - g_val);
      b +=       (1 - a_val) * g_val;
      c +=       a_val * (1 - g_val);
      d +=             a_val * g_val;
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
dispatch::cpu::getORperm(NumericVector& out, const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm) {
  CHECK_VECTOR_SIZE(out, 8); // Check if output vector has correct size

  int ncol = A.ncol();
  int a = 0, b = 0, c = 0, d = 0;

  for (int j = 0; j < ncol - 1; ++j) {
    for (int i = j + 1; i < ncol; ++i) {
      int a_val = A(perm[i], perm[j]);
      int g_val = G(i, j);
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
dispatch::cpu::permuteOR(NumericMatrix& out, const IntegerMatrix& A, const IntegerMatrix& G, int p) {
  CHECK_MATRIX_DIMS(out, p, 8);
  int ncol = A.ncol();
  NumericVector or_tmp(8);
  for (int i = 0; i < p; ++i) {
    IntegerVector perm = sample(ncol, ncol, false) - 1;
    dispatch::cpu::getORperm(or_tmp, A, G, perm);
    for (int j = 0; j < 8; ++j) {
        out(i, j) = or_tmp[j];
    }
  }
}

void
dispatch::cpu::getFDR(List& out, double actual, const NumericVector& permuted) {
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
dispatch::cpu::getG(IntegerMatrix& out, const IntegerVector& Gk) {
  int n = Gk.size();
  CHECK_MATRIX_DIMS(out, n, n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      int gi = Gk[i];
      out(i, i) = gi * gi;
      int value = gi * Gk[j];
      out(i, j) = value;
      out(j, i) = value;
    }
  }
}

void
dispatch::cpu::graflex(NumericVector& out, const IntegerMatrix& A, const IntegerVector& Gk, int p) {
  IntegerMatrix G_tmp(Gk.size(), Gk.size());
  dispatch::cpu::getG(G_tmp, Gk);

  NumericVector actual_tmp(8);
  dispatch::cpu::getOR(actual_tmp, A, G_tmp);

  CHECK_VECTOR_SIZE(out, actual_tmp.length());
  for (int i = 0; i < actual_tmp.length(); ++i) {
      out[i] = actual_tmp[i];
  }

  if (!std::isnan(actual_tmp(4))) {
    NumericMatrix permuted_tmp(p, 8);
    dispatch::cpu::permuteOR(permuted_tmp, A, G_tmp, p);

    List fdr_tmp; 
    dispatch::cpu::getFDR(fdr_tmp, actual_tmp(4), permuted_tmp(_, 4));
    out(6) = as<double>(fdr_tmp["under"]);
    out(7) = as<double>(fdr_tmp["over"]);
  }
}