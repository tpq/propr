#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
void setCudaProfile(bool enable) { (void)enable; }

// [[Rcpp::export]]
DataFrame consumeCudaProfile() {
    return DataFrame::create(
        _["name"] = CharacterVector(0),
        _["ms"] = NumericVector(0),
        _["stringsAsFactors"] = false);
}
