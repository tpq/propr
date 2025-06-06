#include <Rcpp.h>
#include <math.h>
#include "../../../include/kernels/cpu/dispatch/lr2propr.hpp"
#include "../../../include/kernels/cpu/dispatch/backend.hpp"
#include "../../../include/utils.hpp"

using namespace Rcpp;
using namespace propr;


void
dispatch::cpu::lr2vlr(Rcpp::NumericMatrix &lr, NumericMatrix& out) {
  int nfeats = lr.ncol();
  CHECK_MATRIX_DIMS(out, nfeats, nfeats);

  NumericMatrix x_copy = clone(lr);
  NumericMatrix cov_tmp(nfeats, nfeats);
  dispatch::cpu::covRcpp(x_copy, 0, cov_tmp);
  for (int i = 0; i < nfeats; ++i) {
      for (int j = 0; j < nfeats; ++j) {
          out(i, j) = cov_tmp(i, j);
      }
  }

  NumericVector diag(nfeats);
  for(int j = 0; j < nfeats; j++){
    diag[j] = out(j, j);
  }
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      out(i, j) = -2 * out(i, j) + diag[i] + diag[j];
    }
  }
}

void
dispatch::cpu::lr2phi(Rcpp::NumericMatrix &lr, NumericMatrix& out) {
  int nfeats = lr.ncol();
  NumericMatrix x_copy = clone(lr);
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2vlr(x_copy, mat_tmp);
  CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  for (int i = 0; i < mat_tmp.nrow(); ++i) {
      for (int j = 0; j < mat_tmp.ncol(); ++j) {
          out(i, j) = mat_tmp(i, j);
      }
  }

  int nsubjs = lr.nrow();
  for(int i = 0; i < out.ncol(); i++){
    double vari = sum(pow(lr(_, i) - mean(lr(_, i)), 2.0)) / (nsubjs - 1);
    out(_, i) = out(_, i) / vari;
    out(i, i) = 0;
  }
}

void
dispatch::cpu::lr2rho(Rcpp::NumericMatrix &lr, NumericMatrix& out) {
  int nfeats = lr.ncol();
  NumericMatrix x_copy = clone(lr);
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2vlr(x_copy, mat_tmp);
  CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  for (int i = 0; i < mat_tmp.nrow(); ++i) {
      for (int j = 0; j < mat_tmp.ncol(); ++j) {
          out(i, j) = mat_tmp(i, j);
      }
  }
  int nsubjs = lr.nrow();
  NumericVector vars(nfeats);
  for(int i = 0; i < nfeats; i++){
    NumericVector lr_col = lr(_, i);
    vars[i] = sum(pow(lr_col - mean(lr_col), 2.0)) / (nsubjs - 1);
  }
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      if(i == j){
        out(i, j) = 1;
      }else{
        out(i, j) = 1 - out(i, j) / (vars[i] + vars[j]);
      }
    }
  }
}

void
dispatch::cpu::lr2phs(Rcpp::NumericMatrix &lr, NumericMatrix& out) {
  int nfeats = lr.ncol();
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2rho(lr, mat_tmp);
  CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  for (int i = 0; i < mat_tmp.nrow(); ++i) {
      for (int j = 0; j < mat_tmp.ncol(); ++j) {
          out(i, j) = mat_tmp(i, j);
      }
  }
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      if(i == j){
        out(i, j) = 0;
      }else{
        if(out(i, j) == 0){
          out(i, j) = R_PosInf;
        }else{
          out(i, j) = (1 - out(i, j)) / (1 + out(i, j));
        }
      }
    }
  }
}