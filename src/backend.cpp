#include <Rcpp.h>
#include <math.h>
#include <stdint.h>
#include "backend.h"
using namespace Rcpp;

// Calculate weighted mean
// [[Rcpp::export]]
double wtmRcpp(NumericVector x,
               NumericVector w){

  return sum(x * w) / sum(w);
}

// Calculate weighted var
// [[Rcpp::export]]
double wtvRcpp(NumericVector x,
               NumericVector w){

  double xbar = wtmRcpp(x, w);
  return sum(w * pow(x - xbar, 2)) * (sum(w) / (pow(sum(w), 2) - sum(pow(w, 2))));
}

// Function for centering matrix (via: correlateR package)
NumericMatrix centerNumericMatrix(NumericMatrix & X){

  const int32_t m = X.ncol();
  for (int32_t j = 0; j < m; ++j) {
    X(_, j) = X(_, j) - mean(X(_, j));
  }

  return X;
}

// Function for cor (via: correlateR package)
// [[Rcpp::export]]
NumericMatrix corRcpp(NumericMatrix & X){

  const int32_t m = X.ncol();

  // Centering the matrix
  X = centerNumericMatrix(X);

  // Compute 1 over the sample standard deviation
  NumericVector inv_sqrt_ss(m);
  for (int32_t i = 0; i < m; ++i) {
    inv_sqrt_ss(i) = 1 / sqrt(sum(X(_, i) * X(_, i)));
  }

  // Computing the correlation matrix
  NumericMatrix cor(m, m);
  for (int32_t i = 0; i < m; ++i) {
    for (int32_t j = 0; j <= i; ++j) {
      cor(i, j) = sum(X(_,i) * X(_,j)) *
        inv_sqrt_ss(i) * inv_sqrt_ss(j);
      cor(j, i) = cor(i, j);
    }
  }

  return cor;
}

// Function for cov (via: correlateR package)
// [[Rcpp::export]]
NumericMatrix covRcpp(NumericMatrix & X,
                      const int norm_type = 0){

  const int32_t n = X.nrow();
  const int32_t m = X.ncol();
  const int32_t df = n - 1 + norm_type;

  // Centering the matrix
  X = centerNumericMatrix(X);

  // Computing the covariance matrix
  NumericMatrix cov(m, m);
  for (int32_t i = 0; i < m; ++i) {
    for (int32_t j = 0; j <= i; ++j) {
      cov(i, j) = sum(X(_, i) * X(_, j)) / df;
      cov(j, i) = cov(i, j);
    }
  }

  return cov;
}

// Function for vlr
// [[Rcpp::export]]
NumericMatrix vlrRcpp(NumericMatrix & X){

  // Define total number of features
  int32_t nfeats = X.ncol();

  // Take log of "count matrix"
  for(int32_t i = 0; i < X.nrow(); i++){
    for(int32_t j = 0; j < nfeats; j++){
      X(i, j) = log(X(i, j));
    }
  }

  // Calculate variation matrix
  X = covRcpp(X);

  // Find diagonal
  NumericVector diag(nfeats);
  for(int32_t j = 0; j < nfeats; j++){
    diag[j] = X(j, j);
  }

  // Calculate vlr
  for(int32_t i = 0; i < nfeats; i++){
    for(int32_t j = 0; j < nfeats; j++){
      X(i, j) = -2 * X(i, j) + diag[i] + diag[j];
    }
  }

  return X;
}

// Function for clr
// [[Rcpp::export]]
NumericMatrix clrRcpp(NumericMatrix & X){

  // Take log of "count matrix"
  for(int32_t i = 0; i < X.nrow(); i++){
    for(int32_t j = 0; j < X.ncol(); j++){
      X(i, j) = log(X(i, j));
    }

    // Subtract out row mean from logX
    X(i, _) = X(i, _) - mean(X(i, _));
  }

  return X;
}

// Function for alr
// [[Rcpp::export]]
NumericMatrix alrRcpp(NumericMatrix & X,
                      const int32_t ivar = 0){

  if(ivar == 0){

    stop("Select non-zero ivar.");
  }

  // Take log of "count matrix"
  for(int32_t i = 0; i < X.nrow(); i++){
    for(int32_t j = 0; j < X.ncol(); j++){
      X(i, j) = log(X(i, j));
    }

    // Subtract out row ivar from logX
    X(i, _) = X(i, _) - X(i, ivar - 1);
  }

  return X;
}

// Function for matrix symmetry
// [[Rcpp::export]]
NumericMatrix symRcpp(NumericMatrix & X){

  for(int32_t i = 1; i < X.nrow(); i++){
    for(int32_t j = 0; j < i; j++){
      X(j, i) = X(i, j);
    }
  }

  return X;
}

// Function for phi
// [[Rcpp::export]]
NumericMatrix phiRcpp(NumericMatrix X,
                      const bool sym = 1){

  // Make vlr from "count matrix" X (using a copy)
  NumericMatrix counts = clone(X);
  NumericMatrix mat = vlrRcpp(counts);

  // Make clr from "count matrix" X
  NumericMatrix clr = clrRcpp(X);
  int32_t nsubjs = clr.nrow();

  for(int32_t i = 0; i < mat.ncol(); i++){

    // Calculate variance of the i-th clr composition
    double vari = sum(pow(clr(_, i) - mean(clr(_, i)), 2.0)) / (nsubjs - 1);

    // Calculate phi
    mat(_, i) = mat(_, i) / vari;
  }

  if(sym){

    mat = symRcpp(mat);
  }

  return mat;
}

// Function for rho
// [[Rcpp::export]]
NumericMatrix rhoRcpp(NumericMatrix X,
                      NumericMatrix lr,
                      const int32_t ivar = 0){

  // Make vlr from "count matrix" X (using a copy)
  NumericMatrix counts = clone(X);
  NumericMatrix mat = vlrRcpp(counts);
  int32_t nsubjs = X.nrow();
  int32_t nfeats = X.ncol();

  // Calculate variance of the i-th clr composition
  NumericVector vars(nfeats);
  for(int32_t i = 0; i < nfeats; i++){

    vars[i] = sum(pow(lr(_, i) - mean(lr(_, i)), 2.0)) / (nsubjs - 1);
  }

  for(int32_t i = 0; i < nfeats; i++){
    for(int32_t j = 0; j < nfeats; j++){
      if(i == (ivar - 1) || j == (ivar - 1)){

        // Set rho = 0 when ivar is row or column
        if(i == (ivar - 1) && j == (ivar - 1)){
          mat(i, j) = 1;
        }else{
          mat(i, j) = 0;
        }

      }else{

        // Calculate rho
        mat(i, j) = 1 - mat(i, j) / (vars[i] + vars[j]);
      }
    }
  }

  return mat;
}

// Function for pair indexing
// [[Rcpp::export]]
std::vector<int32_t> indexPairs(NumericMatrix & X,
                                const String op = "==",
                                const double ref = 0){

  // Index pairs by operator and reference
  std::vector<int32_t> index;
  int32_t nfeats = X.nrow();
  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){

      if(op == "==" || op == "="){
        if(X(i, j) == ref){
          index.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">"){
        if(X(i, j) > ref){
          index.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">="){
        if(X(i, j) >= ref){
          index.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<"){
        if(X(i, j) < ref){
          index.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<="){
        if(X(i, j) <= ref){
          index.push_back(j * nfeats + i + 1);
        }
      }else if(op == "!="){
        if(X(i, j) != ref){
          index.push_back(j * nfeats + i + 1);
        }

      }else if(op == "all"){

        index.push_back(j * nfeats + i + 1);

      }else{

        stop("Operator not found.");
      }
    }
  }

  return index;
}

// Function for pair indexing
// [[Rcpp::export]]
List indexToCoord(IntegerVector V, const int32_t N){

  std::vector<int32_t> rows;
  std::vector<int32_t> cols;

  for(int32_t i = 0; i < V.length(); i++){

    int32_t J = (V[i] - 1) / N + 1;
    cols.push_back(J);
    int32_t I = (V[i] - 1) % N + 1;
    rows.push_back(I);
  }

  return List::create(
    _["feat1"] = rows,
    _["feat2"] = cols
  );
}

// Function for pair indexing
// [[Rcpp::export]]
IntegerVector coordToIndex(IntegerVector row,
                           IntegerVector col,
                           const int32_t N){

  IntegerVector V = (col - 1) * N + row;
  return V;
}

// Function for Lin's Z and its variance
// [[Rcpp::export]]
NumericMatrix linRcpp(NumericMatrix & rho,
                      NumericMatrix lr){

  int32_t N = lr.nrow();
  NumericMatrix r = corRcpp(lr);

  for(int32_t i = 1; i < rho.nrow(); i++){
    for(int32_t j = 0; j < i; j++){

      // Calculate Z and variance of Z
      double var_ij = (1 - pow(r(i, j), 2)) * pow(rho(i, j), 2) /
      (1 - pow(rho(i, j), 2)) / pow(r(i, j), 2) / (N - 2);
      double z_ij = atanh(rho(i, j));

      // Replace r with Z and var
      r(j, i) = var_ij;
      r(i, j) = z_ij;
    }
  }

  return r;
}

// Function to return lower left triangle
// [[Rcpp::export]]
NumericVector lltRcpp(NumericMatrix & X){

  int32_t nfeats = X.nrow();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      result(counter) = X(i, j);
      counter += 1;
    }
  }

  return result;
}

// Function to return upper right triangle
// [[Rcpp::export]]
NumericVector urtRcpp(NumericMatrix & X){

  int32_t nfeats = X.nrow();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      result(counter) = X(j, i);
      counter += 1;
    }
  }

  return result;
}

// Function to label a half matrix
// [[Rcpp::export]]
List labRcpp(int32_t nfeats){

  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::IntegerVector partner(llt);
  Rcpp::IntegerVector pair(llt);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      partner(counter) = i + 1;
      pair(counter) = j + 1;
      counter += 1;
    }
  }

  return List::create(
    _["Partner"] = partner,
    _["Pair"] = pair
  );
}

// Function to make matrix from half-matrix
// [[Rcpp::export]]
NumericMatrix half2mat(NumericVector X){

  int32_t nfeats = sqrt(2 * X.length() + .25) + .5;
  NumericMatrix A(nfeats, nfeats);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      A(i, j) = X(counter);
      A(j, i) = X(counter);
      counter += 1;
    }
  }

  return A;
}

// Function to recast matrix as feature ratios
// [[Rcpp::export]]
NumericMatrix ratiosRcpp(NumericMatrix & X){

  int32_t nfeats = X.ncol();
  int32_t nsamps = X.nrow();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericMatrix result(nsamps, llt);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){
      result(_, counter) = X(_, i) / X(_, j);
      counter += 1;
    }
  }

  return result;
}
