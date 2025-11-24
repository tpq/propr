#include <Rcpp.h>
#include <math.h>
#include <propr/kernels/cpu/dispatch/lr2propr.hpp>
#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/host_exclusive_profiler.hpp>

using namespace Rcpp;
using namespace propr;


void
dispatch::cpu::lr2vlr(NumericMatrix& out, Rcpp::NumericMatrix &lr) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int nfeats = lr.ncol();
  PROPR_CHECK_MATRIX_DIMS(out, nfeats, nfeats);

  NumericMatrix x_copy = clone(lr);
  NumericMatrix cov_tmp(nfeats, nfeats);
  dispatch::cpu::covRcpp(out, x_copy, 0);
  
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
dispatch::cpu::lr2phi(NumericMatrix& out, Rcpp::NumericMatrix &lr) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int nfeats = lr.ncol();
  NumericMatrix x_copy = clone(lr);
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2vlr(mat_tmp, x_copy);
  PROPR_CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  
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
dispatch::cpu::lr2rho(NumericMatrix& out, Rcpp::NumericMatrix &lr) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int nfeats = lr.ncol();
  NumericMatrix x_copy = clone(lr);
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2vlr(mat_tmp, x_copy);
  PROPR_CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
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
dispatch::cpu::lr2phs(NumericMatrix& out, Rcpp::NumericMatrix &lr) {
  PROPR_PROFILE_HOST_EXCLUSIVE("kernel"); 
  int nfeats = lr.ncol();
  NumericMatrix mat_tmp(nfeats, nfeats);
  dispatch::cpu::lr2rho(mat_tmp, lr);
  PROPR_CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
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