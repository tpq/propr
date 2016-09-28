#include <Rcpp.h>

using namespace Rcpp;

// Function for centering matrix (via: correlateR package)
NumericMatrix centerNumericMatrix(NumericMatrix & X){

  const int m = X.ncol();
  for (int j = 0; j < m; ++j) {
    X(_, j) = X(_, j) - mean(X(_, j));
  }

  return X;
}

// Function for cov (via: correlateR package)
// [[Rcpp::export]]
NumericMatrix covRcpp(NumericMatrix & X,
                      const int norm_type = 0){

  const int n = X.nrow();
  const int m = X.ncol();
  const int df = n - 1 + norm_type;

  // Centering the matrix!
  X = centerNumericMatrix(X);  // Defined in aux_functions

  // Computing the covariance matrix
  NumericMatrix cov(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
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
  int nfeats = X.ncol();

  // Take log of "count matrix"
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < nfeats; j++){
      X(i, j) = log(X(i, j));
    }
  }

  // Calculate variation matrix
  X = covRcpp(X);

  // Find diagonal
  NumericVector diag(nfeats);
  for(int j = 0; j < nfeats; j++){
    diag[j] = X(j, j);
  }

  // Calculate vlr
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      X(i, j) = -2 * X(i, j) + diag[i] + diag[j];
    }
  }

  return X;
}

// Function for clr
// [[Rcpp::export]]
NumericMatrix clrRcpp(NumericMatrix & X){

  // Take log of "count matrix"
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
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
                      const int ivar = 0){

  if(ivar == 0){

    stop("Select non-zero ivar.");
  }

  // Take log of "count matrix"
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
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

  for(int i = 1; i < X.nrow(); i++){
    for(int j = 0; j < i; j++){
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
  int nsubjs = clr.nrow();

  for(int i = 0; i < mat.ncol(); i++){

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
                      const int ivar = 0){


  // Make vlr from "count matrix" X (using a copy)
  NumericMatrix counts = clone(X);
  NumericMatrix mat = vlrRcpp(counts);

  // Create lr matrix container
  int nsubjs = X.nrow();
  int nfeats = X.ncol();
  NumericMatrix lr(nsubjs, nfeats);

  if(ivar == 0){

    // Make clr from "count matrix" X
    lr = clrRcpp(X);

  }else{

    // Make alr from "count matrix" X
    lr = alrRcpp(X, ivar);
  }

  // Calculate variance of the i-th clr composition
  NumericVector vars(nfeats);
  for(int i = 0; i < nfeats; i++){

    vars[i] = sum(pow(lr(_, i) - mean(lr(_, i)), 2.0)) / (nsubjs - 1);
  }

  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      if(i == (ivar - 1) | j == (ivar - 1)){

        // Set rho = 0 when ivar is row or column
        if(i == (ivar - 1) & j == (ivar - 1)){
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
std::vector<int> indexPairs(NumericMatrix & X,
                            const String op = "==",
                            const double ref = 0){

  // Index pairs by operator and reference
  std::vector<int> index;
  int nfeats = X.nrow();
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){

      if(op == "==" | op == "="){
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
List indexToCoord(IntegerVector V, const int N){

  std::vector<int> rows;
  std::vector<int> cols;

  for(int i = 0; i < V.length(); i++){

    int J = (V[i] - 1) / N + 1;
    cols.push_back(J);
    int I = (V[i] - 1) % N + 1;
    rows.push_back(I);
  }

  return List::create(
    _["feat1"] = rows,
    _["feat2"] = cols
  );
}

/*** R
X <- t(data.frame("a" = sample(1:10), "b" = sample(1:10), "c" = sample(1:10),
                  "d" = sample(1:10), "e" = sample(1:10), "f" = sample(1:10)))

if(!all(round(cov(X) - covRcpp(X), 5) == 0)) stop("covRcpp error!")
if(!all(round(propr:::proprVLR(X) - vlrRcpp(X), 5) == 0)) stop("vlrRcpp error!")
if(!all(round(propr:::proprCLR(X) - clrRcpp(X), 5) == 0)) stop("clrRcpp error!")
if(!all(round(propr:::proprALR(X, ivar = 5) - alrRcpp(X, ivar = 5)[, -5], 5) == 0)) stop("alrRcpp error!")
if(!all((propr:::proprSym(X) - symRcpp(X)) == 0)) stop("symRcpp error!")

if(!all(round(propr:::proprPhit(X) - phiRcpp(X), 5) == 0)) stop("phiRcpp error!")
if(!all(round(propr:::proprPerb(X) - rhoRcpp(X), 5) == 0)) stop("rhoRcpp error!")
if(!all(propr:::proprTri(X) - X[indexPairs(X, "all")] == 0)) stop("indexPairs error!")

*/
