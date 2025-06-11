#include <propr/interface/device_selector.hpp>
#include <Rcpp.h>

bool propr::is_gpu_backend() {
    static bool printed = false;
    Rcpp::Environment base = Rcpp::Environment::base_env();
    Rcpp::Function getOption = base["getOption"];
    Rcpp::RObject result = getOption(Rcpp::CharacterVector::create("propr.use_gpu"));    
    if (result.isNULL()) {
        if (!printed) {
            printed = true;
            Rprintf("propr.use_gpu = FALSE (default)\n");
        }
        return false;
    }    
    bool use_gpu = Rcpp::as<bool>(result);
    if (!printed) {
        printed = true;
        Rprintf("propr.use_gpu = %s\n", (use_gpu ? "TRUE" : "FALSE"));
    }
    return use_gpu;
}