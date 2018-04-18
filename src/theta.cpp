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
NumericVector lrm(NumericMatrix & X,
                  NumericMatrix & W,
                  bool weighted = false){

  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;

  if(weighted){
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

  return result;
}

// Calculate lrv with or without weights
// [[Rcpp::export]]
NumericVector lrv(NumericMatrix & Y,
                  NumericMatrix & W,
                  bool weighted = false,
                  double a = NA_REAL){

  // Output a half-matrix
  NumericMatrix X = clone(Y);
  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector result(llt);
  int counter = 0;

  if(!R_IsNA(a)){ // Weighted and non-weighted, alpha-transformed

    // Raise all of X to the a power
    for(int i = 0; i < X.nrow(); i++){
      for(int j = 0; j < X.ncol(); j++){
        X(i, j) = pow(X(i, j), a);
      }
    }

    if(weighted){

      // Sweep out (weighted) column means
      // Calculate sum(W * [i - j]^2)
      // Divide sum(W * [i - j]^2) by (p * a^2)
      Rcpp::NumericVector Wij(nfeats);
      Rcpp::NumericVector Xiscaled(nfeats);
      Rcpp::NumericVector Xjscaled(nfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          Xiscaled = X(_, i) / wtmRcpp(X(_, i), Wij);
          Xjscaled = X(_, j) / wtmRcpp(X(_, j), Wij);
          result(counter) = sum(Wij * pow(Xiscaled - Xjscaled, 2));
          result(counter) = result(counter) /
            (pow(a, 2) * (sum(Wij) - sum(pow(Wij, 2)) / sum(Wij)));
          counter += 1;
        }
      }

    }else{

      // Sweep out column means
      for(int j = 0; j < X.ncol(); j++){
        X(_, j) = X(_, j) / mean(X(_, j));
      }

      // Calculate sum([i - j]^2)
      // Divide sum([i - j]^2) by ((N-1) * a^2)
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          result(counter) = sum(pow(X(_, i) - X(_, j), 2));
          result(counter) = result(counter) /
            (pow(a, 2) * (X.nrow() - 1));
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

// Calculate lrv w.r.t. z
// [[Rcpp::export]]
NumericVector lrz(NumericMatrix & Y,
                  NumericMatrix & W,
                  NumericVector & Z,
                  bool weighted = false,
                  double a = NA_REAL){

  // Output a half-matrix
  NumericMatrix X = clone(Y);
  NumericVector z = clone(Z);
  int nfeats = X.ncol();
  int nsamps = X.nrow();
  NumericVector result(nfeats);
  int counter = 0;

  if(!R_IsNA(a)){ // Weighted and non-weighted, alpha-transformed

    // Raise all of X and z to the a power
    for(int i = 0; i < nsamps; i++){
      z(i) = pow(Z(i), a);
      for(int j = 0; j < nfeats; j++){
        X(i, j) = pow(X(i, j), a);
      }
    }

    if(weighted){

      // Sweep out (weighted) column means
      // Calculate sum(W * [i - j]^2)
      // Divide sum(W * [i - j]^2) by (p * a^2)
      Rcpp::NumericVector Wij(nfeats);
      Rcpp::NumericVector Xiscaled(nfeats);
      Rcpp::NumericVector Xjscaled(nfeats);
      for(int j = 0; j < nfeats; j++){
        Wij = W(_, j);
        Xiscaled = X(_, j) / wtmRcpp(X(_, j), Wij);
        Xjscaled = z / wtmRcpp(z, Wij);
        result(counter) = sum(Wij * pow(Xiscaled - Xjscaled, 2));
        result(counter) = result(counter) /
          (pow(a, 2) * (sum(Wij) - sum(pow(Wij, 2)) / sum(Wij)));
        counter += 1;
      }

    }else{

      // Sweep out column means
      z = z / mean(z);
      for(int j = 0; j < nfeats; j++){
        X(_, j) = X(_, j) / mean(X(_, j));
      }

      // Calculate sum([i - j]^2)
      // Divide sum([i - j]^2) by ((N-1) * a^2)
      for(int j = 0; j < nfeats; j++){
        result(counter) = sum(pow(X(_, j) - z, 2));
        result(counter) = result(counter) /
          (pow(a, 2) * (X.nrow() - 1));
        counter += 1;
      }
    }


  }else{ // Weighted and non-weighted, non-transformed

    if(weighted){

      Rcpp::NumericVector Wij(nfeats);
      for(int j = 0; j < nfeats; j++){
        Wij = W(_, j);
        result(counter) = wtvRcpp(log(X(_, j) / z), Wij);
        counter += 1;
      }

    }else{

      for(int j = 0; j < nfeats; j++){
        result(counter) = var(log(X(_, j) / z));
        counter += 1;
      }
    }
  }

  return result;
}
