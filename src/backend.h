#ifndef BACKEND
#define BACKEND

#include <Rcpp.h>
#include <stdint.h>
#include <math.h>

using namespace Rcpp;

NumericMatrix covRcpp(NumericMatrix & X, const int32_t norm_type);
double wtmRcpp(NumericVector x, NumericVector w);
double wtvRcpp(NumericVector x, NumericVector w);

#endif
