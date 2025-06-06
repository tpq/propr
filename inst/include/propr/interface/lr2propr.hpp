#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::NumericMatrix lr2vlr(Rcpp::NumericMatrix lr);
    Rcpp::NumericMatrix lr2phi(Rcpp::NumericMatrix lr);
    Rcpp::NumericMatrix lr2rho(Rcpp::NumericMatrix lr);
    Rcpp::NumericMatrix lr2phs(Rcpp::NumericMatrix lr);
}