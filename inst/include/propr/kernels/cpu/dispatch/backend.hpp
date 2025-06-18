#pragma once
#include <Rcpp.h>

namespace propr {
    namespace dispatch {
        namespace cpu {
            
            __forceinline__
            void wtmRcpp(double& out, const Rcpp::NumericVector &x, Rcpp::NumericVector &w);
            
            __forceinline__
            void wtvRcpp(double& out, const Rcpp::NumericVector &x, Rcpp::NumericVector &w);
            
            __forceinline__
            void corRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void covRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X, int norm_type = 0);
            
            __forceinline__
            void vlrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void clrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void alrRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X, int ivar = 0);
            
            __forceinline__
            void symRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void phiRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &X, bool sym = true);
            
            __forceinline__
            void rhoRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix &X, Rcpp::NumericMatrix& lr, int ivar = 0);
            
            __forceinline__
            void indexPairs(std::vector<int>& out, Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0);
            
            __forceinline__
            void indexToCoord(Rcpp::List& out, Rcpp::IntegerVector V, int N);
            
            __forceinline__
            void coordToIndex(Rcpp::IntegerVector& out, Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N);
            
            __forceinline__
            void linRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix & rho, Rcpp::NumericMatrix lr);
            
            __forceinline__
            void lltRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void urtRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void labRcpp(Rcpp::List & out, int nfeats);
            
            __forceinline__
            void half2mat(Rcpp::NumericMatrix& out, Rcpp::NumericVector X);
            
            __forceinline__
            void vector2mat(Rcpp::NumericMatrix& out, Rcpp::NumericVector X, Rcpp::IntegerVector i, Rcpp::IntegerVector j, int nfeats);
            
            __forceinline__
            void ratiosRcpp(Rcpp::NumericMatrix & out, Rcpp::NumericMatrix & X);
            
            __forceinline__
            void results2matRcpp(Rcpp::NumericMatrix & out, Rcpp::DataFrame& results, int n, double diagonal = 0.0);
        }
    }
}