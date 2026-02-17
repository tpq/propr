#include <Rcpp.h>

using namespace Rcpp;

void setCudaProfile(bool enable) { (void)enable; }

DataFrame consumeCudaProfile() {
    return DataFrame::create(
        _["name"] = CharacterVector(0),
        _["ms"] = NumericVector(0),
        _["stringsAsFactors"] = false);
}
