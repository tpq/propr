#pragma once
#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            void dof_global    (Rcpp::NumericVector& out, const Rcpp::NumericMatrix & W, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void dof_population(Rcpp::NumericVector& out, const Rcpp::NumericMatrix & W, propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}