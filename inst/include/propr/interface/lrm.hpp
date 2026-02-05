#pragma once

#include <Rcpp.h>

namespace propr {
	Rcpp::NumericVector lrm(Rcpp::NumericMatrix &Y,
	                        Rcpp::NumericMatrix &W,
	                        bool weighted = false,
	                        double a = NA_REAL,
	                        Rcpp::NumericMatrix Yfull = Rcpp::NumericMatrix(1, 1),
	                        Rcpp::NumericMatrix Wfull = Rcpp::NumericMatrix(1, 1),
	                        bool use_gpu = false);
}
