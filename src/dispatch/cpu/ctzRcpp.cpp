#include <Rcpp.h>

#include <propr/utils/rcpp_checks.h>
#include <propr/kernels/cpu/dispatch/ctzRcpp.hpp>


using namespace Rcpp;
using namespace propr;

void
dispatch::cpu::ctzRcpp(NumericVector& out, NumericMatrix & X){
    int nfeats = X.ncol();
    int nsubjs = X.nrow();
    int llt = nfeats * (nfeats - 1) / 2;
    CHECK_VECTOR_SIZE(out, llt);
    NumericVector zeroes(nfeats);
    for(int i = 0; i < nfeats; i++){
        for(int j = 0; j < nsubjs; j++){
            if(X(j, i) == 0){
                zeroes(i) += 1;
            }
        }
    }
    int counter = 0;
    for(int i = 1; i < nfeats; i++){
        for(int j = 0; j < i; j++){
            out(counter) = zeroes(i) + zeroes(j);
            counter += 1;
        }
    }
}