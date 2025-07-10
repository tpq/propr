#include <Rcpp.h>
#include <numeric>

#include <propr/utils/rcpp_checks.h>
#include <propr/kernels/cpu/dispatch/omega.hpp>

using namespace Rcpp;

void
propr::dispatch::cpu::dof_global(NumericVector& out, NumericMatrix & W) {
	int nfeats = W.ncol();
	int llt = nfeats * (nfeats - 1) / 2;
	CHECK_VECTOR_SIZE(out, llt);
	NumericVector Wij(W.nrow());
	int counter = 0;
	double n    = 0;
	for(int i = 1; i < nfeats; i++){
		for(int j = 0; j < nfeats; j++){
			Wij = W(_, i) * W(_, j);
			n = sum(Wij);
			out(counter) = n - sum(pow(Wij, 2)) / n;
			counter += 1;
		}
	}
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