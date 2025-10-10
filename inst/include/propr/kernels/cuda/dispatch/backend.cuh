#pragma once
#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            void wtmRcpp(double& out, const Rcpp::NumericVector &x, const Rcpp::NumericVector &w, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void wtvRcpp(double& out, const Rcpp::NumericVector &x, const Rcpp::NumericVector &w, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void corRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void covRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, int norm_type = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void vlrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void clrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void alrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, int ivar = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void symRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void phiRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &X, bool sym = true, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void rhoRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix &X,const Rcpp::NumericMatrix& lr, int ivar = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void indexPairs(std::vector<int>& out, const Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void indexToCoord(Rcpp::List& out, const Rcpp::IntegerVector V, int N, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void coordToIndex(Rcpp::IntegerVector& out, const Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void linRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & rho, const Rcpp::NumericMatrix& lr, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void lltRcpp(Rcpp::NumericVector& out,const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void urtRcpp(Rcpp::NumericVector& out,const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void labRcpp(Rcpp::List & out, int nfeats, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void half2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void vector2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X, const Rcpp::IntegerVector& i, const Rcpp::IntegerVector& j, int nfeats, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void ratiosRcpp(Rcpp::NumericMatrix & out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            void results2matRcpp(Rcpp::NumericMatrix & out, const Rcpp::DataFrame& results, int n, double diagonal = 0.0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}