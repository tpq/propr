#include <propr/interface/device_selector.hpp>
#include <Rcpp.h>


bool is_gpu_backend() {
    Rcpp::Environment base = Rcpp::Environment::base_env();
    Rcpp::Function getOption = base["getOption"];
    Rcpp::RObject  result = getOption(Rcpp::CharacterVector::create("propr.use.gpu"));
    if (result.isNULL())  return false;
    return Rcpp::as<bool>(result);
}