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
	int llt = nfeats * (nfeats - 1) / 2;
	PROPR_CHECK_VECTOR_SIZE(out, llt)
	int counter = 0;

	if (!R_IsNA(a)) {
		if (Yfull.nrow() == NumericMatrix(1, 1).nrow() &&
			Yfull.ncol() == NumericMatrix(1, 1).ncol()) {
			stop("User must provide valid Yfull argument for alpha-transformation.");
		}

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
		if (weighted) {
			if (Wfull.nrow() == NumericMatrix(1, 1).nrow() &&
				Wfull.ncol() == NumericMatrix(1, 1).ncol()) {
				stop("User must provide valid Wfull argument for weighted alpha-transformation.");
			}

			Rcpp::NumericVector Wij(X.nrow());
			Rcpp::NumericVector Wfullij(Wfull.nrow());
			Rcpp::NumericVector Xiscaled(X.nrow());
			Rcpp::NumericVector Xjscaled(X.nrow());
			Rcpp::NumericVector Xz(X.nrow());
			Rcpp::NumericVector Xfullz(Xfull_copy.nrow());

            double mean_Xfull_i, mean_Xfull_j, Mz_sum_val;

			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					Wij = W(_, i) * W(_, j);
					Wfullij = Wfull(_, i) * Wfull(_, j);

                    dispatch::cpu::wtmRcpp(mean_Xfull_i, Xfull_copy(_, i), Wfullij);
                    dispatch::cpu::wtmRcpp(mean_Xfull_j, Xfull_copy(_, j), Wfullij);
					Xiscaled = X(_, i) / mean_Xfull_i;
					Xjscaled = X(_, j) / mean_Xfull_j;

					Xz = X(_, i) - X(_, j);
					Xfullz = Xfull_copy(_, i) - Xfull_copy(_, j);

                    dispatch::cpu::wtmRcpp(Mz_sum_val, Xiscaled - Xjscaled, Wij);
					double Mz = Mz_sum_val;

					double Cz = sum(Wij * Xz) / sum(Wij) + (sum(Wfullij * Xfullz) - sum(Wij * Xz)) / (sum(Wfullij) - sum(Wij));
					out(counter) = (Cz / 2 + Mz) / a;
					counter += 1;
				}
			}
		} else {
			Rcpp::NumericVector Xz(X.nrow());
			Rcpp::NumericVector Xfullz(Yfull.nrow());
			double N1 = X.nrow();
			double NT = Yfull.nrow();
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					Xz = X(_, i) - X(_, j);
					Xfullz = Yfull(_, i) - Yfull(_, j);

					double Mz = mean(X(_, i) / mean(Yfull(_, i)) - X(_, j) / mean(Yfull(_, j)));
					double Cz = sum(Xz) / N1 + (sum(Xfullz) - sum(Xz)) / (NT - N1);
					out(counter) = (Cz / 2 + Mz) / a;
					counter += 1;
				}
			}
		}
	}
	else { // Weighted and non-weighted, non-transformed
		if (weighted) {
			Rcpp::NumericVector Wij(X.nrow());
            double lrm_val;
			for (int i = 1; i < nfeats; i++) {
				for (int j = 0; j < i; j++) {
					Wij = W(_, i) * W(_, j);
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
