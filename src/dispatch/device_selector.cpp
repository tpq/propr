#include <propr/interface/device_selector.hpp>
#include <Rcpp.h>

bool propr::is_gpu_backend() {
    static bool printed = false;
    Rcpp::Environment base = Rcpp::Environment::base_env();
    Rcpp::Function getOption = base["getOption"];
    Rcpp::RObject  result = getOption(Rcpp::CharacterVector::create("propr.use_gpu"));
    if (!printed){
        printed = true;
        bool use_gpu = Rcpp::as<bool>(result);
        Rprintf("propr.use_gpu = %s\n", (use_gpu ? "TRUE" : "FALSE"));
    }
    if (result.isNULL())  return false;
    return Rcpp::as<bool>(result);
}