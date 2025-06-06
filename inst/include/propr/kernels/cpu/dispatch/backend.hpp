#pragma once
#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void wtmRcpp(const Rcpp::NumericVector &x, Rcpp::NumericVector &w,double& out);
            void wtvRcpp(const Rcpp::NumericVector &x, Rcpp::NumericVector &w, double& out);
            void corRcpp(Rcpp::NumericMatrix & X,  Rcpp::NumericMatrix& out);
            void covRcpp(Rcpp::NumericMatrix & X, int norm_type = 0,  Rcpp::NumericMatrix& out);
            void vlrRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out);
            void clrRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out);
            void alrRcpp(Rcpp::NumericMatrix & X, int ivar = 0,  Rcpp::NumericMatrix& out);
            void symRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericMatrix& out);
            void phiRcpp(Rcpp::NumericMatrix &X, bool sym = true,  Rcpp::NumericMatrix& out);
            void rhoRcpp(Rcpp::NumericMatrix &X, Rcpp::NumericMatrix& lr, int ivar = 0,  Rcpp::NumericMatrix& out);
            void indexPairs(Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0, std::vector<int>& out);
            void indexToCoord(Rcpp::IntegerVector V, int N, Rcpp::List& out);
            void coordToIndex(Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, Rcpp::IntegerVector& out);
            void linRcpp(Rcpp::NumericMatrix & rho, Rcpp::NumericMatrix lr, Rcpp::NumericMatrix& out);
            void lltRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericVector& out);
            void urtRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericVector& out);
            void labRcpp(int nfeats, Rcpp::List & out);
            void half2mat(Rcpp::NumericVector X, Rcpp::NumericMatrix& out);
            void vector2mat(Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats, Rcpp::NumericMatrix& out);
            void ratiosRcpp(Rcpp::NumericMatrix & X, Rcpp::NumericMatrix & out);
            void results2matRcpp(Rcpp::DataFrame& results, int n, double diagonal = 0.0, Rcpp::NumericMatrix & out );
        }
    }
}