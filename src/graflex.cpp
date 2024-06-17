// #include <RcppArmadillo.h>
#include <Rcpp.h>
#include <numeric>

using namespace Rcpp;

// // Function to extract the lower triangle of a matrix
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::mat getLowerTriangle(const arma::mat& mat) {
//     return arma::trimatl(mat);
// }

// Function to set the seed for reproducibility
void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// Function to extract the lower triangle of a square IntegerMatrix
IntegerVector get_square_integer_matrix_triangle(const IntegerMatrix& mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  int n = ncol * (ncol - 1) / 2;
  IntegerVector triangle(n);

  if (nrow != ncol) {
    stop("Input matrix must be square");
  }

  int k = 0;
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < i; ++j) {
      triangle[k++] = mat(i, j);
    }
  }

  return triangle;
}

// Function to shuffle a square IntegerMatrix by column and row
// [[Rcpp::export]]
IntegerMatrix shuffle_square_integer_matrix(const IntegerMatrix& mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  IntegerMatrix shuffled(nrow, ncol);

  if (nrow != ncol) {
    stop("Input matrix must be square");
  }

  // shuffle by index
  IntegerVector index = sample(ncol, ncol, false) - 1;

  // fill in the shuffled matrix
  // NOTE that the diagonal stays 0
  for (int i = 0; i < ncol; ++i) {
    for (int j = 0; j < i; ++j) {
      shuffled(i, j) = mat(index[i], index[j]);
      shuffled(j, i) = mat(index[j], index[i]);
    }
  }

  return shuffled;
}

// Function to calculate the odds ratio table
// [[Rcpp::export]]
NumericMatrix binTab(IntegerVector& A, IntegerVector& G) {
  NumericMatrix mat(2, 2);
  int n = A.size();

  mat(1, 1) = sum(A * G);              // in A and G
  mat(1, 0) = sum(A * (1 - G));        // in A but not G
  mat(0, 0) = sum((1 - A) * (1 - G));  // not in A and not in G
  mat(0, 1) = sum((1 - A) * G);        // not in A but in G

  return mat;
}

// [[Rcpp::export]]
NumericMatrix getOR(IntegerVector& A, IntegerVector& G) {
  NumericMatrix tab = binTab(A, G);

  double odds_ratio = (tab(0, 0) * tab(1, 1)) / (tab(0, 1) * tab(1, 0));

  NumericMatrix res(1, 6);
  res(0, 0) = tab(0, 0);
  res(0, 1) = tab(0, 1);
  res(0, 2) = tab(1, 0);
  res(0, 3) = tab(1, 1);
  res(0, 4) = odds_ratio;
  res(0, 5) = log(odds_ratio);

  return res;
}

// Function to calculate the odds ratio and other relevant info for each permutation
// [[Rcpp::export]]
NumericMatrix permuteOR(IntegerMatrix A, IntegerMatrix G, int p = 100, Nullable<int> seed = R_NilValue) {
  NumericMatrix or_table(p, 6);
  int nrow = A.nrow();
  int ncol = A.ncol();

  if (nrow != ncol) {
    stop("Input matrix must be square");
  }
  if (nrow != G.nrow() || ncol != G.ncol()) {
    stop("Input matrices must have the same dimensions");
  }

  if (seed.isNotNull()) set_seed(Rcpp::as<int>(seed));

  IntegerVector Gstar = get_square_integer_matrix_triangle(G);

  // calculate odds ratio for each permutation
  for (int i = 0; i < p; ++i) {
    IntegerMatrix Ashuffled = shuffle_square_integer_matrix(A);
    IntegerVector Astar = get_square_integer_matrix_triangle(Ashuffled);
    or_table(i, _) = getOR(Astar, Gstar);
  }

  return(or_table);
}

