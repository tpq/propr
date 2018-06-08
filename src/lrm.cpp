#include <Rcpp.h>
#include <math.h>
#include "backend.h"
using namespace Rcpp;

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
