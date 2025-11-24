#include <Rcpp.h>
#include <cmath>
#include <propr/kernels/cpu/dispatch/backend.hpp>
#include <propr/kernels/cpu/dispatch/lrv.hpp>
#include <propr/utils/rcpp_checks.h>
#include <propr/utils/host_profiler.hpp>

using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::lrv(Rcpp::NumericVector& out,
				  Rcpp::NumericMatrix &Y,
				  Rcpp::NumericMatrix &W,
				  bool weighted,
				  double a,
				  Rcpp::NumericMatrix Yfull,
				  Rcpp::NumericMatrix Wfull) {
	NumericMatrix X = clone(Y);
	int nfeats        = X.ncol();
	int fullfeats     = Yfull.ncol();
	int llt           = nfeats * (nfeats - 1) / 2;
	PROPR_CHECK_VECTOR_SIZE(out, llt);
	int counter = 0;

	if (!R_IsNA(a)) {

		if (Yfull.nrow() == NumericMatrix(1, 1).nrow() &&
			Yfull.ncol() == NumericMatrix(1, 1).ncol()) {
			stop("User must provide valid Yfull argument for alpha-transformation.");
		}

		if (weighted) {
			PROPR_PROFILE_HOST("kernel"); 

			for (int i = 0; i < X.nrow(); i++) {
				for (int j = 0; j < X.ncol(); j++) {
					X(i, j) = pow(X(i, j), a);
				}
			}

			NumericMatrix Xfull_copy = clone(Yfull);
			int full_X_rows = Xfull_copy.nrow();
			int full_X_cols = Xfull_copy.ncol();
			for (int i = 0; i < full_X_rows; i++) {
				for (int j = 0; j < full_X_cols; j++) {
					Xfull_copy(i, j) = pow(Xfull_copy(i, j), a);
				}
			}

			if (Wfull.nrow() == NumericMatrix(1, 1).nrow() &&
				Wfull.ncol() == NumericMatrix(1, 1).ncol()) {
				stop("User must provide valid Wfull argument for weighted alpha-transformation.");
			}

            double mean_X_i, mean_X_j, mean_Xfull_i, mean_Xfull_j;
			// Mean-center the within-group values as a fraction of the across-group means
			// Calculate sum(W * [i - j]^2)
			// Divide sum(W * [i - j]^2) by (p * a^2)
			Rcpp::NumericVector Wij(nfeats);
     		Rcpp::NumericVector Wfullij(fullfeats);
     		Rcpp::NumericVector Xiscaled(nfeats);
     		Rcpp::NumericVector Xjscaled(nfeats);

			for(int i = 1; i < nfeats; i++){
				for(int j = 0; j < i; j++){
					Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
					Wfullij = 2 * Wfull(_, i) * Wfull(_, j) / (Wfull(_, i) + Wfull(_, j));

					dispatch::cpu::wtmRcpp(mean_X_i, X(_, i), Wij);
					dispatch::cpu::wtmRcpp(mean_X_j, X(_, j), Wij);
					dispatch::cpu::wtmRcpp(mean_Xfull_i,Xfull_copy(_, i), Wfullij);
					dispatch::cpu::wtmRcpp(mean_Xfull_j,Xfull_copy(_, j), Wfullij);

					Xiscaled = (X(_, i) - mean_X_i) / mean_Xfull_i;
					Xjscaled = (X(_, j) - mean_X_j) / mean_Xfull_j;
					out(counter) = sum(Wij * pow(Xiscaled - Xjscaled, 2)) / (pow(a, 2) * (sum(Wij) - sum(pow(Wij, 2)) / sum(Wij)));
					counter += 1;
				}
			}
		} else {
			PROPR_PROFILE_HOST("kernel"); 
			for (int i = 0; i < X.nrow(); i++) {
				for (int j = 0; j < X.ncol(); j++) {
					X(i, j) = pow(X(i, j), a);
				}
			}

			NumericMatrix Xfull_copy = clone(Yfull);
			int full_X_rows = Xfull_copy.nrow();
			int full_X_cols = Xfull_copy.ncol();
			for (int i = 0; i < full_X_rows; i++) {
				for (int j = 0; j < full_X_cols; j++) {
					Xfull_copy(i, j) = pow(Xfull_copy(i, j), a);
				}
			}
			Rcpp::NumericVector Xiscaled(nfeats);
			Rcpp::NumericVector Xjscaled(nfeats);
			double N1 = X.nrow();

            double mean_Xfull_i, mean_Xfull_j;

			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
                    mean_Xfull_i = mean(Xfull_copy(_, i));
                    mean_Xfull_j = mean(Xfull_copy(_, j));

					Xiscaled = (X(_, i) - mean(X(_, i))) / mean_Xfull_i;
					Xjscaled = (X(_, j) - mean(X(_, j))) / mean_Xfull_j;

					out(counter) = sum(pow(Xiscaled - Xjscaled, 2)) /(pow(a, 2) * (N1 - 1));
					counter += 1;
				}
			}
		}
	} else { // Weighted and non-weighted, non-transformed
		if (weighted) {
			PROPR_PROFILE_HOST("kernel"); 
			NumericVector Wij(X.nrow());
            double wtv_val;
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					Wij = 2 * W(_, i) * W(_, j) / (W(_, i) + W(_, j));
					dispatch::cpu::wtvRcpp(wtv_val, log(X(_, i) / X(_, j)), Wij);
					out(counter) = wtv_val;
					counter += 1;
				}
			}
		} else {
			PROPR_PROFILE_HOST("kernel"); 
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					out(counter) = var(log(X(_, i) / X(_, j)));
					counter += 1;
				}
			}
		}
	}
}