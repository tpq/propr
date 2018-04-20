#include <Rcpp.h>

using namespace Rcpp;

// Function to count joint zero frequency
// [[Rcpp::export]]
NumericVector ctzRcpp(NumericMatrix & X){

  int nfeats = X.ncol();
  int nsubjs = X.nrow();
  int llt = nfeats * (nfeats - 1) / 2;

  // Count zero frequency per feature
  Rcpp::NumericVector zeroes(nfeats);
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nsubjs; j++){
      if(X(j, i) == 0){
        zeroes(i) += 1;
      }
    }
  }

  // Count joint zero frequency
  Rcpp::NumericVector result(llt);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      result(counter) = zeroes(i) + zeroes(j);
      counter += 1;
    }
  }

  return result;
}

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

// Calculate lrm with or without weights
// [[Rcpp::export]]
NumericVector lrm(NumericMatrix & Y,
                  NumericMatrix & W,
                  bool weighted = false,
                  double a = NA_REAL,
                  NumericMatrix Yfull = NumericMatrix(1, 1),
                  NumericMatrix Wfull = NumericMatrix(1, 1)){

  // Output a half-matrix
  NumericMatrix X = clone(Y);
  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector result(llt);
  int counter = 0;

  if(!R_IsNA(a)){ // Weighted and non-weighted, alpha-transformed

    //////////////////////////////////////////////////////
    // Check for valid Yfull argument
    if(Yfull.nrow() == NumericMatrix(1, 1).nrow() &&
       Yfull.ncol() == NumericMatrix(1, 1).ncol()){

      stop("User must provide valid Yfull argument for alpha-transformation.");
    }

    // Raise all of X to the a power
    for(int i = 0; i < X.nrow(); i++){
      for(int j = 0; j < X.ncol(); j++){
        X(i, j) = pow(X(i, j), a);
      }
    }

    // Raise all of Xfull to the a power
    NumericMatrix Xfull = clone(Yfull);
    int fullfeats = Xfull.ncol();
    for(int i = 0; i < Xfull.nrow(); i++){
      for(int j = 0; j < Xfull.ncol(); j++){
        Xfull(i, j) = pow(Xfull(i, j), a);
      }
    }
    //////////////////////////////////////////////////////

    if(weighted){

      //////////////////////////////////////////////////////
      // Check for valid Wfull argument
      if(Wfull.nrow() == NumericMatrix(1, 1).nrow() &&
         Wfull.ncol() == NumericMatrix(1, 1).ncol()){

        stop("User must provide valid Wfull argument for weighted alpha-transformation.");
      }
      //////////////////////////////////////////////////////

      // Calculate alpha-transformed mean using the across-group means
      Rcpp::NumericVector Wij(nfeats);
      Rcpp::NumericVector Wfullij(fullfeats);
      Rcpp::NumericVector Xiscaled(nfeats);
      Rcpp::NumericVector Xjscaled(nfeats);
      Rcpp::NumericVector Xz(nfeats);
      Rcpp::NumericVector Xfullz(fullfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          Wfullij = Wfull(_, i) * Wfull(_, j);
          Xiscaled = X(_, i) / wtmRcpp(Xfull(_, i), Wfullij);
          Xjscaled = X(_, j) / wtmRcpp(Xfull(_, j), Wfullij);
          Xz = X(_, i) - X(_, j);
          Xfullz = Xfull(_, i) - Xfull(_, j);
          double Mz = sum(Wij * (Xiscaled - Xjscaled) / sum(Wij));
          double Cz = sum(Wij * Xz) / sum(Wij) +
            (sum(Wfullij * Xfullz) - sum(Wij * Xz)) /
              (sum(Wfullij) - sum(Wij));
          result(counter) = (Cz/2 + Mz) / a;
          counter += 1;
        }
      }

    }else{

      // Calculate alpha-transformed mean using the across-group means
      Rcpp::NumericVector Xz(nfeats);
      Rcpp::NumericVector Xfullz(fullfeats);
      double N1 = X.nrow();
      double NT = Xfull.nrow();
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Xz = X(_, i) - X(_, j);
          Xfullz = Xfull(_, i) - Xfull(_, j);
          double Mz = mean(X(_, i) / mean(Xfull(_, i)) - mean(X(_, j) / mean(Xfull(_, j))));
          double Cz = sum(Xz) / N1 +
            (sum(Xfullz) - sum(Xz)) /
              (NT - N1);
          result(counter) = (Cz/2 + Mz) / a;
          counter += 1;
        }
      }
    }

  }else{ // Weighted and non-weighted, non-transformed

    if(weighted){

      Rcpp::NumericVector Wij(nfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          result(counter) = wtmRcpp(log(X(_, i) / X(_, j)), Wij);
          counter += 1;
        }
      }

    }else{

      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          result(counter) = mean(log(X(_, i) / X(_, j)));
          counter += 1;
        }
      }
    }
  }

  return result;
}

// Calculate lrv with or without weights
// [[Rcpp::export]]
NumericVector lrv(NumericMatrix & Y,
                  NumericMatrix & W,
                  bool weighted = false,
                  double a = NA_REAL,
                  NumericMatrix Yfull = NumericMatrix(1, 1),
                  NumericMatrix Wfull = NumericMatrix(1, 1)){

  // Output a half-matrix
  NumericMatrix X = clone(Y);
  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector result(llt);
  int counter = 0;

  if(!R_IsNA(a)){ // Weighted and non-weighted, alpha-transformed

    //////////////////////////////////////////////////////
    // Check for valid Yfull argument
    if(Yfull.nrow() == NumericMatrix(1, 1).nrow() &&
       Yfull.ncol() == NumericMatrix(1, 1).ncol()){

      stop("User must provide valid Yfull argument for alpha-transformation.");
    }

    // Raise all of X to the a power
    for(int i = 0; i < X.nrow(); i++){
      for(int j = 0; j < X.ncol(); j++){
        X(i, j) = pow(X(i, j), a);
      }
    }

    // Raise all of Xfull to the a power
    NumericMatrix Xfull = clone(Yfull);
    int fullfeats = Xfull.ncol();
    for(int i = 0; i < Xfull.nrow(); i++){
      for(int j = 0; j < Xfull.ncol(); j++){
        Xfull(i, j) = pow(Xfull(i, j), a);
      }
    }
    //////////////////////////////////////////////////////

    if(weighted){

      //////////////////////////////////////////////////////
      // Check for valid Wfull argument
      if(Wfull.nrow() == NumericMatrix(1, 1).nrow() &&
         Wfull.ncol() == NumericMatrix(1, 1).ncol()){

        stop("User must provide valid Wfull argument for weighted alpha-transformation.");
      }
      //////////////////////////////////////////////////////

      // Mean-center the within-group values as a fraction of the across-group means
      // Calculate sum(W * [i - j]^2)
      // Divide sum(W * [i - j]^2) by (p * a^2)
      Rcpp::NumericVector Wij(nfeats);
      Rcpp::NumericVector Wfullij(fullfeats);
      Rcpp::NumericVector Xiscaled(nfeats);
      Rcpp::NumericVector Xjscaled(nfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          Wfullij = Wfull(_, i) * Wfull(_, j);
          Xiscaled = (X(_, i) - wtmRcpp(X(_, i), Wij)) / wtmRcpp(Xfull(_, i), Wfullij);
          Xjscaled = (X(_, j) - wtmRcpp(X(_, j), Wij)) / wtmRcpp(Xfull(_, j), Wfullij);
          result(counter) = sum(Wij * pow(Xiscaled - Xjscaled, 2)) /
            (pow(a, 2) * (sum(Wij) - sum(pow(Wij, 2)) / sum(Wij)));
          counter += 1;
        }
      }

    }else{

      // Mean-center the within-group values as a fraction of the across-group means
      // Calculate sum([i - j]^2)
      // Divide sum([i - j]^2) by ((N-1) * a^2)
      Rcpp::NumericVector Xiscaled(nfeats);
      Rcpp::NumericVector Xjscaled(nfeats);
      double N1 = X.nrow();
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Xiscaled = (X(_, i) - mean(X(_, i))) / mean(Xfull(_, i));
          Xjscaled = (X(_, j) - mean(X(_, j))) / mean(Xfull(_, j));
          result(counter) = sum(pow(Xiscaled - Xjscaled, 2)) /
            (pow(a, 2) * (N1 - 1));
          counter += 1;
        }
      }
    }

  }else{ // Weighted and non-weighted, non-transformed

    if(weighted){

      Rcpp::NumericVector Wij(nfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          result(counter) = wtvRcpp(log(X(_, i) / X(_, j)), Wij);
          counter += 1;
        }
      }

    }else{

      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          result(counter) = var(log(X(_, i) / X(_, j)));
          counter += 1;
        }
      }
    }
  }

  return result;
}

// Calculate lrv weight modifier
// [[Rcpp::export]]
NumericVector omega(NumericMatrix & X,
                    NumericMatrix & W){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;
  double n = 0;

  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      Wij = W(_, i) * W(_, j);
      n = sum(Wij);
      result(counter) = n - sum(pow(Wij, 2)) / n;
      counter += 1;
    }
  }

  return result;
}
