#pragma once
#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void wtmRcpp(double& out, const Rcpp::NumericVector &x, Rcpp::NumericVector &w);
            void wtvRcpp(double& out, const Rcpp::NumericVector &x, Rcpp::NumericVector &w);
            void corRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            void covRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X, int norm_type = 0);
            void vlrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            void clrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            void alrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X, int ivar = 0);
            void symRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            void phiRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &X, bool sym = true);
            void rhoRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &X, Rcpp::NumericMatrix& lr, int ivar = 0);
            void indexPairs(std::vector<int>& out, Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0);
            void indexToCoord(Rcpp::List& out, Rcpp::IntegerVector V, int N);
            void coordToIndex(Rcpp::IntegerVector& out, Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N);
            void linRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & rho, Rcpp::NumericMatrix lr);
            void lltRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix & X);
            void urtRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix & X);
            void labRcpp(Rcpp::List & out, int nfeats);
            void half2mat(Rcpp::NumericMatrix& out, Rcpp::NumericVector X);
            void vector2mat(Rcpp::NumericMatrix& out, Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats);
            void ratiosRcpp(Rcpp::NumericMatrix & out, Rcpp::NumericMatrix & X);
            void results2matRcpp(Rcpp::NumericMatrix & out, Rcpp::DataFrame& results, int n, double diagonal = 0.0);
        }
    }
}