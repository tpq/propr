#include <Rcpp.h>
#include <numeric>

#include <propr/utils/rcpp_checks.h>
#include <propr/kernels/cpu/dispatch/omega.hpp>

using namespace Rcpp;


NumericMatrix pad_matrix2(const NumericMatrix& mat,
                         int padTop, int padBottom,
                         int padLeft, int padRight) {
  int old_nrow = mat.nrow();
  int old_ncol = mat.ncol();
  int new_nrow = old_nrow + padTop + padBottom;
  int new_ncol = old_ncol + padLeft + padRight;
  
  // create output matrix initialized to zero
  NumericMatrix out(new_nrow, new_ncol);
  
  // copy original into the right offset
  for (int i = 0; i < old_nrow; ++i) {
    for (int j = 0; j < old_ncol; ++j) {
      out(i + padTop, j + padLeft) = mat(i, j);
    }
  }
  
  return out;
}

void
propr::dispatch::cpu::dof_global(NumericVector& out, NumericMatrix & W) {

  int nfeats = W.ncol();
  int llt = nfeats * (nfeats - 1) / 2;
  Rcpp::NumericVector result(llt);
  Rcpp::NumericVector Wij(nfeats);
  int counter = 0;
  double n = 0;

  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){
      Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
      n = sum(Wij);
      out(counter) = n - sum(pow(Wij, 2)) / n;
      std::cout << out(counter) << " ";
      counter += 1;
    } std::cout << std::endl;
  }
}

void
propr::dispatch::cpu::dof_population(NumericVector& out, const NumericMatrix & W) {
    int nfeats = W.ncol();
    int llt    = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    int counter = 0;

    NumericVector Wij(nfeats);
    for(int i = 1; i < nfeats; i++){
      for(int j = 0; j < i; j++){
        Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
        out(counter) = sum(Wij);
        counter += 1;
      }
    }
}
