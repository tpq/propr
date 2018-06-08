#ifndef BACKEND
#define BACKEND

#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

NumericMatrix covRcpp(NumericMatrix & X, const int norm_type);
double wtmRcpp(NumericVector x, NumericVector w);
double wtvRcpp(NumericVector x, NumericVector w);

#endif
