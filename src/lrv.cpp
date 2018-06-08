#include <Rcpp.h>
#include <math.h>
#include "backend.h"
using namespace Rcpp;

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
