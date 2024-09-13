#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;


// Function to extract the triangle of a square and symmetric IntegerMatrix
// [[Rcpp::export]]
IntegerVector get_triangle(const IntegerMatrix& mat) {
  int ncol = mat.ncol();
  int n = ncol * (ncol - 1) / 2;

  IntegerVector triangle(n);

  int k = 0;
  for (int j = 0; j < ncol; ++j) {
      for (int i = j+1; i < ncol; ++i) {
      triangle[k++] = mat(i, j);
    }
  }

  return triangle;
}

// Function to get the triangle of a matrix based on the given indeces
// [[Rcpp::export]]
IntegerVector get_triangle_from_index(const IntegerMatrix& mat, const IntegerVector& index) {
  int ncol = mat.ncol();
  int n = ncol * (ncol - 1) / 2;

  IntegerVector triangle(n);

  int k = 0;
  for (int j = 0; j < ncol; ++j) {
    for (int i = j+1; i < ncol; ++i) {
      triangle[k++] = mat(index[i], index[j]);
    }
  }

  return triangle;
}

// Function to calculate the contingency table
// [[Rcpp::export]]
NumericVector getOR(const IntegerVector& A, const IntegerVector& G) {
  int n = A.size();

  // calculate the contingency table
  int a = 0, b = 0, c = 0, d = 0;
  for (int i = 0; i < n; ++i) {
    if (A[i] == 0) {
      if (G[i] == 0) ++a;  // not in A and not in G
      else ++b;            // not in A but in G
    } else {
      if (G[i] == 0) ++c;  // in A but not in G
      else ++d;            // in A and in G
    }
  }
  
  // calculate the odds ratio
  double odds_ratio = static_cast<double>(a * d) / (b * c);

  return NumericVector::create(
    a, b, c, d, odds_ratio, std::log(odds_ratio), R_NaN, R_NaN
  );
}

// Function to calculate the odds ratio and other relevant info for each permutation
// [[Rcpp::export]]
NumericMatrix permuteOR(const IntegerMatrix& A, const IntegerVector& Gstar, int p = 100) {
  int n = A.ncol();
  NumericMatrix or_table(p, 8);

  // calculate odds ratio for each permutation
  for (int i = 0; i < p; ++i) {
    IntegerVector idx = sample(n, n, false) - 1;
    IntegerVector Astar = get_triangle_from_index(A, idx);
    or_table(i, _) = getOR(Astar, Gstar);
  }

  return(or_table);
}

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

// Function to calculate the odds ratio and FDR, given the adjacency matrix A and the knowledge graph G
// [[Rcpp::export]]
NumericVector graflex(const IntegerMatrix& A, const IntegerMatrix& G, int p = 100) {

  // get the actual odds ratio
  IntegerVector Astar = get_triangle(A);
  IntegerVector Gstar = get_triangle(G);
  NumericVector actual = getOR(Astar, Gstar);

  // get distribution of odds ratios on permuted data
  NumericMatrix permuted = permuteOR(A, Gstar, p);

  // calculate the FDR
  List fdr = getFDR(actual(4), permuted(_, 4));
  actual(6) = as<double>(fdr["under"]);
  actual(7) = as<double>(fdr["over"]);

  return actual;
}
