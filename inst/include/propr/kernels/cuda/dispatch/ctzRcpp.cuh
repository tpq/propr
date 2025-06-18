#pragma once
#include <Rcpp.h>
#include <propr/context.h>


namespace propr {
    namespace dispatch {
        namespace cuda {
            __forceinline__ 
            void ctzRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix & X, propr::propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}