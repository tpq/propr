#include <Rcpp.h>
#include <numeric>
using namespace Rcpp;

// Function to extract the lower triangle of a square and symmetric IntegerMatrix
// [[Rcpp::export]]
IntegerVector get_lower_triangle(IntegerMatrix& mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int n = ncol * (ncol - 1) / 2;

  IntegerVector triangle(n);

  int k = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < i; ++j) {
      triangle[k++] = mat(i, j);
    }
  }

  return triangle;
}

// Function to shuffle a square and symmetric IntegerMatrix, and get the lower triangle
// [[Rcpp::export]]
IntegerVector shuffle_and_get_lower_triangle(IntegerMatrix& mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int n = ncol * (ncol - 1) / 2;

  IntegerVector shuffled_triangle(n);
  IntegerVector index = sample(ncol, ncol, false) - 1;

  int k = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < i; ++j) {
      shuffled_triangle[k++] = mat(index[i], index[j]);
    }
  }

  return shuffled_triangle;
}

// Function to calculate the contingency table
// [[Rcpp::export]]
NumericMatrix binTab(IntegerVector& A, IntegerVector& G) {
  NumericMatrix tab(2, 2);
  int n = A.size();

  tab(0, 0) = sum((1 - A) * (1 - G));  // not in A and not in G
  tab(0, 1) = sum((1 - A) * G);        // not in A but in G
  tab(1, 0) = sum(A * (1 - G));        // in A but not G
  tab(1, 1) = sum(A * G);              // in A and G

  return tab;
}

// Function to calculate the odds ratio, and parse all information into a matrix
// [[Rcpp::export]]
NumericMatrix getOR(IntegerVector& A, IntegerVector& G) {

  NumericMatrix tab = binTab(A, G);
  double odds_ratio = (tab(0, 0) * tab(1, 1)) / (tab(0, 1) * tab(1, 0));

  NumericMatrix or_table(1, 8);
  or_table(0, 0) = tab(0, 0);
  or_table(0, 1) = tab(0, 1);
  or_table(0, 2) = tab(1, 0);
  or_table(0, 3) = tab(1, 1);
  or_table(0, 4) = odds_ratio;
  or_table(0, 5) = log(odds_ratio);

  return or_table;
}

// Function to calculate the odds ratio and other relevant info for each permutation
// [[Rcpp::export]]
NumericMatrix permuteOR(IntegerMatrix& A, IntegerVector& Gstar, int p = 100) {

  NumericMatrix or_table(p, 8);

  // calculate odds ratio for each permutation
  for (int i = 0; i < p; ++i) {
    IntegerVector Astar = shuffle_and_get_lower_triangle(A);
    or_table(i, _) = getOR(Astar, Gstar);
  }

  return(or_table);
}

// [[Rcpp::export]]
float getFDR_over(float actual, NumericVector permuted) {
  float n = permuted.size();
  float fdr = 0.0;
  for (int i = 0; i < n; ++i) {
    if (permuted[i] >= actual) {
      fdr += 1.0;
    }
  }
  fdr /= n;
  return fdr;
}

// [[Rcpp::export]]
float getFDR_under(float actual, NumericVector permuted) {
  float n = permuted.size();
  float fdr = 0.0;
  for (int i = 0; i < n; ++i) {
    if (permuted[i] <= actual) {
      fdr += 1.0;
    }
  }
  fdr /= n;
  return fdr;
}

// Function to calculate the odds ratio and FDR, given the adjacency matrix A and the knowledge graph G
// [[Rcpp::export]]
NumericMatrix graflex(IntegerMatrix& A, IntegerMatrix& G, int p = 100) {

  // get the actual odds ratio
  IntegerVector Astar = get_lower_triangle(A);
  IntegerVector Gstar = get_lower_triangle(G);
  NumericMatrix actual = getOR(Astar, Gstar);

  // get distribution of odds ratios on permuted data
  NumericMatrix permuted = permuteOR(A, Gstar, p);

  // calculate the FDR
  actual(0, 6) = getFDR_under(actual(0, 4), permuted(_, 4));
  actual(0, 7) = getFDR_over(actual(0, 4), permuted(_, 4));

  return actual;
}
