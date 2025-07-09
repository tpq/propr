#include <Rcpp.h>
#include <math.h>
#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/utils.hpp>

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::wtmRcpp(double& out, const NumericVector& x, NumericVector& w){
  out =  sum(x * w) / sum(w); 
}

void
dispatch::cpu::wtvRcpp(double& out, const NumericVector& x, NumericVector& w) {
  double xbar;
  wtmRcpp(xbar, x, w);
  out = sum(w * pow(x - xbar, 2)) * (sum(w) / (pow(sum(w), 2) - sum(pow(w, 2))));
}

void
centerNumericMatrix(NumericMatrix& out, NumericMatrix & in_X){
  CHECK_MATRIX_DIMS(out, in_X.nrow(), in_X.ncol());
  for (int j = 0; j < in_X.ncol(); ++j) {
    for (int i = 0; i < in_X.nrow(); ++i) {
      out(i, j) = in_X(i, j);
    }
  }
  const int m = out.ncol();
  for (int j = 0; j < m; ++j) {
    out(_, j) = out(_, j) - mean(out(_, j));
  }
}

void
dispatch::cpu::corRcpp(NumericMatrix& out, NumericMatrix & X) {
  const int m = X.ncol();
  CHECK_MATRIX_DIMS(out, m, m);

  NumericMatrix X_centered;
  centerNumericMatrix(X, X_centered);

  NumericVector inv_sqrt_ss(m);
  for (int i = 0; i < m; ++i) {
    inv_sqrt_ss(i) = 1 / sqrt(sum(X_centered(_, i) * X_centered(_, i)));
  }
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      out(i, j) = sum(X_centered(_,i) * X_centered(_,j)) * inv_sqrt_ss(i) * inv_sqrt_ss(j);
      out(j, i) = out(i, j);
    }
  }
}

void
dispatch::cpu::covRcpp(NumericMatrix& out, NumericMatrix & X,const int norm_type) {

  const int n = X.nrow();
  const int m = X.ncol();
  CHECK_MATRIX_DIMS(out, m, m);
  const int df = n - 1 + norm_type;

  NumericMatrix X_centered_temp(n, m);
  centerNumericMatrix(X, X_centered_temp);

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      out(i, j) = sum(X_centered_temp(_, i) * X_centered_temp(_, j)) / df;
      out(j, i) = out(i, j);
    }
  }
}

void
dispatch::cpu::vlrRcpp(NumericMatrix& out, NumericMatrix & X){

  int nfeats = X.ncol();
  CHECK_MATRIX_DIMS(out, nfeats, nfeats);

  NumericMatrix X_log(X.nrow(), X.ncol());
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < nfeats; j++){
      X_log(i, j) = log(X(i, j));
    }
  }

  NumericMatrix cov_matrix(nfeats, nfeats);
  dispatch::cpu::covRcpp(cov_matrix, X_log, 0);

  NumericVector diag(nfeats);
  for(int j = 0; j < nfeats; j++){
    diag[j] = cov_matrix(j, j);
  }

  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      out(i, j) = -2 * cov_matrix(i, j) + diag[i] + diag[j];
    }
  }
}

void
dispatch::cpu::clrRcpp(NumericMatrix& out, NumericMatrix & X) {
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      out(i, j) = log(X(i, j));
    }
    out(i, _) = out(i, _) - mean(out(i, _));
  }
}

void
dispatch::cpu::alrRcpp(NumericMatrix& out, NumericMatrix & X,const int ivar){
  if(ivar == 0) {
    stop("Select non-zero ivar.");
  }
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      out(i, j) = log(X(i, j));
    }
    out(i, _) = out(i, _) - out(i, ivar - 1);
  }
}

void
dispatch::cpu::symRcpp(NumericMatrix& out, NumericMatrix & X) {
  CHECK_MATRIX_DIMS(out, X.nrow(), X.ncol());
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < X.ncol(); j++){
      out(i,j) = X(i,j);
    }
  }
  for(int i = 1; i < out.nrow(); i++){
    for(int j = 0; j < i; j++){
      out(j, i) = out(i, j);
    }
  }
}

void
dispatch::cpu::phiRcpp(NumericMatrix& out, NumericMatrix& X, bool sym) {
  NumericMatrix mat_tmp(X.nrow(), X.ncol());
  dispatch::cpu::vlrRcpp(mat_tmp, X);

  CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  for (int i = 0; i < mat_tmp.nrow(); ++i) {
      for (int j = 0; j < mat_tmp.ncol(); ++j) {
          out(i, j) = mat_tmp(i, j);
      }
  }

  NumericMatrix clr_matrix(X.nrow(), X.ncol());
  dispatch::cpu::clrRcpp(clr_matrix, X);

  int nsubjs = X.nrow();
  for(int i = 0; i < out.ncol(); i++){
    double vari = sum(pow(clr_matrix(_, i) - mean(clr_matrix(_, i)), 2.0)) / (nsubjs - 1);
    out(_, i) = out(_, i) / vari;
    out(i, i) = 0;
  }
  if(sym) {
      NumericMatrix sym_mat_tmp(out.nrow(), out.ncol());
      dispatch::cpu::symRcpp(sym_mat_tmp, out);
      for (int i = 0; i < out.nrow(); ++i) {
          for (int j = 0; j < out.ncol(); ++j) {
              out(i, j) = sym_mat_tmp(i, j);
          }
      }
  }
}

void
dispatch::cpu::rhoRcpp(NumericMatrix& out, NumericMatrix& X, NumericMatrix& lr,int ivar){

  NumericMatrix mat_tmp(X.nrow(), X.ncol());
  dispatch::cpu::vlrRcpp(mat_tmp, X);

  CHECK_MATRIX_DIMS(out, mat_tmp.nrow(), mat_tmp.ncol());
  for (int i = 0; i < mat_tmp.nrow(); ++i) {
      for (int j = 0; j < mat_tmp.ncol(); ++j) {
          out(i, j) = mat_tmp(i, j);
      }
  }

  int nsubjs = X.nrow();
  int nfeats = X.ncol();

  NumericVector vars(nfeats);
  for(int i = 0; i < nfeats; i++){
    vars[i] = sum(pow(lr(_, i) - mean(lr(_, i)), 2.0)) / (nsubjs - 1);
  }
  for(int i = 0; i < nfeats; i++){
    for(int j = 0; j < nfeats; j++){
      if(i == (ivar - 1) || j == (ivar - 1)){
        if(i == (ivar - 1) && j == (ivar - 1)){
          out(i, j) = 1;
        }else{
          out(i, j) = 0;
        }
      }else{
        out(i, j) = 1 - out(i, j) / (vars[i] + vars[j]);
      }
    }
  }
}

void
dispatch::cpu::indexPairs(std::vector<int>& out, NumericMatrix & X, const String op,const double ref) {
  out.clear();
  int nfeats = X.nrow();
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      if(op == "==" || op == "="){
        if(X(i, j) == ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">"){
        if(X(i, j) > ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == ">="){
        if(X(i, j) >= ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<"){
        if(X(i, j) < ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "<="){
        if(X(i, j) <= ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "!="){
        if(X(i, j) != ref){
          out.push_back(j * nfeats + i + 1);
        }
      }else if(op == "all"){
        out.push_back(j * nfeats + i + 1);
      }else{
        stop("Operator not found.");
      }
    }
  }
}

void
dispatch::cpu::indexToCoord(List& out, IntegerVector V, const int N){
  std::vector<int> rows;
  std::vector<int> cols;
  for(int i = 0; i < V.length(); i++){
    int J = (V[i] - 1) / N + 1;
    cols.push_back(J);
    int I = (V[i] - 1) % N + 1;
    rows.push_back(I);
  }
  out["feat1"] = Rcpp::wrap(rows);
  out["feat2"] = Rcpp::wrap(cols);
}

void
dispatch::cpu::coordToIndex(IntegerVector& out, IntegerVector row,IntegerVector col, const int N){
  CHECK_VECTOR_SIZE(out, row.length());
  if (row.length() != col.length()){
    stop("Input row and col vectors must have the same length.");
  }
  for (int i = 0; i < row.length(); ++i) {
      out[i] = (col[i] - 1) * N + row[i];
  }
}

void
dispatch::cpu::linRcpp(NumericMatrix& out, NumericMatrix & rho, NumericMatrix lr){
  int N_samples = lr.nrow();
  CHECK_MATRIX_DIMS(out, rho.nrow(), rho.ncol());

  NumericMatrix r_tmp(lr.nrow(), lr.ncol());
  dispatch::cpu::corRcpp(r_tmp, lr);

  for(int i = 0; i < rho.nrow(); i++){
    for(int j = 0; j < rho.ncol(); j++){
      out(i,j) = r_tmp(i,j);
    }
  }

  for(int i = 1; i < rho.nrow(); i++){
    for(int j = 0; j < i; j++){
      // Calculate Z and variance of Z
      double var_ij = (1 - pow(r_tmp(i, j), 2)) * pow(rho(i, j), 2) /
      (1 - pow(rho(i, j), 2)) / pow(r_tmp(i, j), 2) / (N_samples - 2);
      double z_ij = atanh(rho(i, j));
      out(j, i) = var_ij;
      out(i, j) = z_ij;
    }
  }
}

void
dispatch::cpu::lltRcpp(NumericVector& out, NumericMatrix & X){
  int nfeats = X.nrow();
  int llt = nfeats * (nfeats - 1) / 2;
  CHECK_VECTOR_SIZE(out, llt);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      out(counter) = X(i, j);
      counter += 1;
    }
  }
}

void
dispatch::cpu::urtRcpp(NumericVector& out, NumericMatrix & X){
  int nfeats = X.nrow();
  int llt = nfeats * (nfeats - 1) / 2;
  CHECK_VECTOR_SIZE(out, llt);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      out(counter) = X(j, i);
      counter += 1;
    }
  }
}

void
dispatch::cpu::labRcpp(List& out, int nfeats){
  int llt = nfeats * (nfeats - 1) / 2;
  out["Partner"] = Rcpp::IntegerVector(llt);
  out["Pair"] = Rcpp::IntegerVector(llt);

  Rcpp::IntegerVector partner_ref = out["Partner"];
  Rcpp::IntegerVector pair_ref = out["Pair"];

  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      partner_ref(counter) = i + 1;
      pair_ref(counter) = j + 1;
      counter += 1;
    }
  }
}

void
dispatch::cpu::half2mat(NumericMatrix& out, NumericVector X){
  int nfeats = sqrt(2 * X.length() + .25) + .5;
  CHECK_MATRIX_DIMS(out, nfeats, nfeats);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      out(i, j) = X(counter);
      out(j, i) = X(counter);
      counter += 1;
    }
  }
}

void
dispatch::cpu::vector2mat(NumericMatrix& out, NumericVector X, IntegerVector i, IntegerVector j, int nfeats) {
  int nX = X.length();
  int ni = i.length();
  int nj = j.length();
  if (ni != nj){
    stop("i and j must be the same length.");
  }
  if (ni != nX){
    stop("i, j, and X must be the same length.");
  }
  CHECK_MATRIX_DIMS(out, nfeats, nfeats);
  std::memset( REAL(out), 0, sizeof(double) * out.size() );
  for (int counter = 0; counter < ni; counter++){
    out(i[counter]-1, j[counter]-1) = X[counter];
    out(j[counter]-1, i[counter]-1) = X[counter];
  }
}

void dispatch::cpu::ratiosRcpp(NumericMatrix& out, NumericMatrix & X){

  int nfeats = X.ncol();
  int nsamps = X.nrow();
  int llt = nfeats * (nfeats - 1) / 2;
  CHECK_MATRIX_DIMS(out, nsamps, llt);
  int counter = 0;
  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      out(_, counter) = X(_, i) / X(_, j);
      counter += 1;
    }
  }
}

// Function to recast results data frame as gene-gene matrix
void
dispatch::cpu::results2matRcpp(Rcpp::NumericMatrix& out, DataFrame& results, int n, double diagonal){
  CHECK_MATRIX_DIMS(out, n, n);
  int npairs = results.nrows();
  for (int i = 0; i < npairs; i++){
    int row = int(results(i, 0)) - 1;
    int col = int(results(i, 1)) - 1;
    out(row, col) = results(i, 2);
    out(col, row) = results(i, 2);
  }
  for (int i = 0; i < n; i++){
    out(i, i) = diagonal;
  }
}