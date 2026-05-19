#pragma once

#include <Rcpp.h>

namespace propr {
    Rcpp::List genewiseConnectivityRcpp(
        const Rcpp::IntegerVector& partner,
        const Rcpp::IntegerVector& pair,
        const Rcpp::NumericVector& theta,
        const Rcpp::NumericVector& fdr,
        int num_genes,
        double pairwise_fdr = 0.05,
        Rcpp::String backend = "auto");
}
