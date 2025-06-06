#pragma once
#include <Rcpp.h>
#include <cuda_runtime.h>
#include "context.h"

namespace propr {
    namespace dispatch {
        namespace cuda {
            void dof_global(const Rcpp::NumericMatrix & W, Rcpp::NumericVector& out, propr_context context);
            void dof_population(const Rcpp::NumericMatrix & W, Rcpp::NumericVector& out, propr_context context);
        }
    }
}