#include <Rcpp.h>

using namespace Rcpp;

// Function for Box-Cox psuedo-vlr
// [[Rcpp::export]]
NumericVector boxRcpp(NumericMatrix & X,
                      const double a){

  // Raise all of X to the a power
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      X(i, j) = pow(X(i, j), a);
    }
  }

  // Sweep out column means
  for(int j = 0; j < X.ncol(); j++){
    X(_, j) = X(_, j) / mean(X(_, j));
  }

  // Output a half-matrix
  int nfeats = X.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector result(llt);
  int counter = 0;

  // Calculate sum([i - j]^2)
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      result(counter) = sum(pow(X(_, i) - X(_, j), 2));
      counter += 1;
    }
  }

  result = result / (pow(a, 2) * (X.nrow() - 1));
  return result;
}

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

// Calculate lrm with weights
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

// Calculate lrv with weights
// [[Rcpp::export]]
NumericVector lrv(NumericMatrix & X,
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

  return result;
}

// Calculate lrv modifier
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


// [[Rcpp::export]]
NumericVector testRcpp(NumericMatrix & Y,
                       NumericMatrix & W,
                       bool weighted = false,
                       const double a = NA_REAL){

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
      for(int j = 0; j < X.ncol(); j++){
        X(_, j) = X(_, j) / wtmRcpp(X(_, j), W(_, j));
      }

      // Calculate sum(W * [i - j]^2)
      // Divide sum(W * [i - j]^2) by (p * a^2)
      Rcpp::NumericVector Wij(nfeats);
      for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
          Wij = W(_, i) * W(_, j);
          result(counter) = sum(Wij * pow(X(_, i) - X(_, j), 2));
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
