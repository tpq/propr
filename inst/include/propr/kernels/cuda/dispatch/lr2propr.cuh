#pragma once

#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            void lr2vlr(Rcpp::NumericMatrix& lr, Rcpp::NumericMatrix &out, propr::propr_context context);
            void lr2phi(Rcpp::NumericMatrix& lr, Rcpp::NumericMatrix &out, propr::propr_context context);

            void lr2rho(Rcpp::NumericMatrix &lr, Rcpp::NumericMatrix &out, propr_context context);
            void lr2phs(Rcpp::NumericMatrix& lr, Rcpp::NumericMatrix &out, propr::propr_context context);
        }
    }
}
