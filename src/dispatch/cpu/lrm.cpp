#include <Rcpp.h>
#include <math.h>
#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/kernels/cpu/dispatch/lrm.hpp>
#include <propr/utils/rcpp_checks.h>

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::lrm(NumericVector& out,
				   NumericMatrix &Y, NumericMatrix &W,
				   bool weighted,double a,
				   NumericMatrix Yfull,NumericMatrix Wfull)
{
	NumericMatrix X = clone(Y);
	int nfeats      = X.ncol();
	int fullfeats   = Yfull.ncol();
	int llt = nfeats * (nfeats - 1) / 2;
	PROPR_CHECK_VECTOR_SIZE(out, llt)
	int counter = 0;

	if (!R_IsNA(a)) { // Weighted and non-weighted, alpha-transformed
		//////////////////////////////////////////////////////
    	// Check for valid Yfull argument
		if (Yfull.nrow() == NumericMatrix(1, 1).nrow() &&
			Yfull.ncol() == NumericMatrix(1, 1).ncol()) {
			stop("User must provide valid Yfull argument for alpha-transformation.");
		}

		// Raise all of X to the a power
		for (int i = 0; i < X.nrow(); i++) {
			for (int j = 0; j < X.ncol(); j++) {
				X(i, j) = pow(X(i, j), a);
			}
		}

		NumericMatrix Xfull = clone(Yfull);
		int full_X_rows = Xfull.nrow();
		int full_X_cols = Xfull.ncol();

		// Raise all of Xfull to the a power
		for (int i = 0; i < full_X_rows; i++) {
			for (int j = 0; j < full_X_cols; j++) {
				Xfull(i, j) = pow(Xfull(i, j), a);
			}
		}
		if (weighted) {
			if (Wfull.nrow() == NumericMatrix(1, 1).nrow() &&
				Wfull.ncol() == NumericMatrix(1, 1).ncol()) {
				stop("User must provide valid Wfull argument for weighted alpha-transformation.");
			}

			// Calculate alpha-transformed mean using the across-group means
			Rcpp::NumericVector Wij(nfeats);
      		Rcpp::NumericVector Wfullij(fullfeats);
      		Rcpp::NumericVector Xiscaled(nfeats);
      		Rcpp::NumericVector Xjscaled(nfeats);
      		Rcpp::NumericVector Xz(nfeats);
      		Rcpp::NumericVector Xfullz(fullfeats);

            double mean_Xfull_i, mean_Xfull_j, Mz_sum_val;

			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
					Wfullij = 2 * Wfull(_, i) * Wfull(_, j) / (Wfull(_, i) + Wfull(_, j));

                    dispatch::cpu::wtmRcpp(mean_Xfull_i, Xfull(_, i), Wfullij);
                    dispatch::cpu::wtmRcpp(mean_Xfull_j, Xfull(_, j), Wfullij);

					Xiscaled = X(_, i) / mean_Xfull_i;
					Xjscaled = X(_, j) / mean_Xfull_j;

					Xz = X(_, i) - X(_, j);
					Xfullz = Xfull(_, i) - Xfull(_, j);
					double Mz = sum(Wij * (Xiscaled - Xjscaled) / sum(Wij));
					double Cz = sum(Wij * Xz) / sum(Wij) + (sum(Wfullij * Xfullz) - sum(Wij * Xz)) / (sum(Wfullij) - sum(Wij));
					out(counter) = (Cz / 2 + Mz) / a;
					counter += 1;
				}
			}
		} else {
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
					double Cz = sum(Xz) / N1 + (sum(Xfullz) - sum(Xz)) /(NT - N1);
					out(counter) = (Cz / 2 + Mz) / a;
					counter += 1;
				}
			}
		}
	}
	else { // Weighted and non-weighted, non-transformed
		if (weighted) {
			Rcpp::NumericVector Wij(nfeats);
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					double lrm_val;
					Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
                    dispatch::cpu::wtmRcpp(lrm_val, log(X(_, i) / X(_, j)), Wij);
					out(counter) = lrm_val;
					counter += 1;
				}
			}
		} else {
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					out(counter) = mean(log(X(_, i) / X(_, j)));
					counter += 1;
				}
			}
		}
	}
}
