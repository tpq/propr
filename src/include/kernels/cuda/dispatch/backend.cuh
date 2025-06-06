#pragma once
#include <Rcpp.h>
#include "../../../context.h"

namespace propr {
    namespace dispatch {
        namespace cuda {
            void wtmRcpp(const Rcpp::NumericVector &x, const Rcpp::NumericVector &w,double& out, propr_context context);
            void wtvRcpp(const Rcpp::NumericVector &x, const Rcpp::NumericVector &w, double& out, propr_context context);
            void corRcpp(const Rcpp::NumericMatrix & X,  Rcpp::NumericMatrix& out, propr_context context);
            void covRcpp(const Rcpp::NumericMatrix & X, int norm_type = 0,  Rcpp::NumericMatrix& out, propr_context context);
            void vlrRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out, propr_context context);
            void clrRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out, propr_context context);
            void alrRcpp(const Rcpp::NumericMatrix & X, int ivar = 0,  Rcpp::NumericMatrix& out, propr_context context);
            void symRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out, propr_context context);
            void phiRcpp(const Rcpp::NumericMatrix &X, bool sym = true,  Rcpp::NumericMatrix& out, propr_context context);
            void rhoRcpp(const Rcpp::NumericMatrix &X,const Rcpp::NumericMatrix& lr, int ivar = 0,  Rcpp::NumericMatrix& out, propr_context context);
            void indexPairs(const Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0, std::vector<int>& out, propr_context context);
            void indexToCoord(const Rcpp::IntegerVector V, int N, Rcpp::List& out, propr_context context);
            void coordToIndex(const Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, Rcpp::IntegerVector& out, propr_context context);
            void linRcpp(const Rcpp::NumericMatrix & rho, const Rcpp::NumericMatrix& lr, Rcpp::NumericMatrix& out, propr_context context);
            void lltRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericVector& out, propr_context context);
            void urtRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericVector& out, propr_context context);
            void labRcpp(int nfeats, Rcpp::List & out, propr_context context);
            void half2mat(const Rcpp::NumericVector& X,Rcpp::NumericMatrix& out, propr_context context);
            void vector2mat(const Rcpp::NumericVector& X, const Rcpp::IntegerVector& i, const Rcpp::IntegerVector& j, int nfeats, Rcpp::NumericMatrix& out, propr_context context);
            void ratiosRcpp(const Rcpp::NumericMatrix & X, Rcpp::NumericMatrix & out, propr_context context);
            void results2matRcpp(const Rcpp::DataFrame& results, int n, double diagonal = 0.0, Rcpp::NumericMatrix & out, propr_context context);
        }
    }
}