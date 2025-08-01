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
	CHECK_VECTOR_SIZE(out, llt);
	NumericVector Wij(W.nrow());

	int counter = 0;
	double n    = 0;
  std::cout << "[CPU]: ";
	for(int i = 1; i < nfeats; i++){
		for(int j = 0; j < i; j++){
			Wij = 2 * W(_, i) * W(_, j) / ( W(_, i) + W(_, j));
			n = sum(Wij);
			out(counter) = n - sum(pow(Wij, 2)) / n;
      std::cout << out(counter) << " ";
			counter += 1;
		}
	}
  std::cout << std::endl;
}

void
propr::dispatch::cpu::dof_population(NumericVector& out, const NumericMatrix & W) {
    int nfeats = W.ncol();
    int llt    = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    int counter = 0;
    for (int i = 1; i < nfeats; ++i) {
        for (int j = 0; j < i; ++j) {
            Rcpp::NumericVector col_i = W(_, i);
            Rcpp::NumericVector col_j = W(_, j);
            double s = std::inner_product(
                col_i.begin(), col_i.end(),
                col_j.begin(),
                0.0
            );
            out[counter++] = s;
        }
    }
}

// [CPU]: 8785.39 
//        2703.29 2435.75 
//        894.549 862.709 705.633 
// [CPU]: 6910.41 12157.3 6096.84 2401.04 2006.21 2293.52 
// [CPU]: 15856.7 14952.2 8605.14 3322.39 2893.77 3021.38 