#pragma once
#include <Rcpp.h>
#include <propr/context.h>


namespace propr {
    namespace dispatch {
        namespace cuda {
            void ctzRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericVector& out, propr::propr_context context);
        }
    }
}