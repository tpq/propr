#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// Optimized function to calculate the contingency table and the odds ratio
// [[Rcpp::export]]
NumericVector getOR(const IntegerMatrix& A, const IntegerMatrix& G) {
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

  return NumericVector::create(a, b, c, d, odds_ratio, log_odds_ratio, R_NaN, R_NaN);
}

// Optimized function to calculate the contingency table and the odds ratio, given a permuted index vector
// [[Rcpp::export]]
NumericVector getORperm(const IntegerMatrix& A, const IntegerMatrix& G, const IntegerVector& perm) {
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

  return NumericVector::create(a, b, c, d, odds_ratio, log_odds_ratio, R_NaN, R_NaN);
}

// Function to calculate the odds ratio and other relevant info for each permutation
// [[Rcpp::export]]
NumericMatrix permuteOR(const IntegerMatrix& A, const IntegerMatrix& G, int p = 100) {
  int ncol = A.ncol();
  NumericMatrix or_table(p, 8);

  // calculate the odds ratio for each permutation
  for (int i = 0; i < p; ++i) {
    IntegerVector perm = sample(ncol, ncol, false) - 1; 
    or_table(i, _) = getORperm(A, G, perm);
    // TODO should I downsample the pairs (up to a maximum) to be checked?
    // So in this case, we would check how likely we get by chance an OR from the downsampled
    // permuted data that is higher/lower than the OR on the downsampled empirical data
  }

  return(or_table);
}

// Function to calculate the FDR, given the actual odds ratio and the permuted odds ratios
// [[Rcpp::export]]
List getFDR(double actual, const NumericVector& permuted) {
  int n = permuted.size();
  int count_over = 0;
  int count_under = 0;

  // Count values above and below the actual value
  for (int i = 0; i < n; ++i) {
    double current = permuted[i];
    if (current >= actual) ++count_over;
    if (current <= actual) ++count_under;
  }

  // Calculate FDR for both "over" and "under"
  double fdr_over = static_cast<double>(count_over) / n;
  double fdr_under = static_cast<double>(count_under) / n;

  // Return both FDR values as a named list
  return List::create(
    Named("over") = fdr_over,
    Named("under") = fdr_under
  );
}

// Function to calculate the G matrix from the Gk vector
// [[Rcpp::export]]
IntegerMatrix getG(const IntegerVector& Gk) {
  int n = Gk.size();
  IntegerMatrix G(n, n);
  
  for (int i = 0; i < n; ++i) {
    int gi = Gk[i];
    G(i, i) = gi * gi;
    for (int j = 0; j < i; ++j) {
      int value = gi * Gk[j];
      G(i, j) = value;
      G(j, i) = value;
    }
  }
  
  return G;
}

// Function to calculate the odds ratio and FDR, given the adjacency matrix A and the knowledge graph G
// [[Rcpp::export]]
NumericVector graflex(const IntegerMatrix& A, const IntegerVector& Gk, int p = 100) {

  // Calculate Gk
  IntegerMatrix G = getG(Gk);

  // get the actual odds ratio
  NumericVector actual = getOR(A, G);

  // get distribution of odds ratios on permuted data
  NumericMatrix permuted = permuteOR(A, G, p);

  // calculate the FDR
  List fdr = getFDR(actual(4), permuted(_, 4));
  actual(6) = as<double>(fdr["under"]);
  actual(7) = as<double>(fdr["over"]);

  return actual;
}
