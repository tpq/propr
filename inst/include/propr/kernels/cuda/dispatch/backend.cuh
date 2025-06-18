#pragma once
#include <Rcpp.h>
#include <propr/context.h>

namespace propr {
    namespace dispatch {
        namespace cuda {
            
            __forceinline__
            void wtmRcpp(double& out, const Rcpp::NumericVector &x, const Rcpp::NumericVector &w, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void wtvRcpp(double& out, const Rcpp::NumericVector &x, const Rcpp::NumericVector &w, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void corRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void covRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, int norm_type = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void vlrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void clrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void alrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, int ivar = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void symRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void phiRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix &X, bool sym = true, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void rhoRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix &X,const Rcpp::NumericMatrix& lr, int ivar = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void indexPairs(std::vector<int>& out, const Rcpp::NumericMatrix & X,Rcpp::String op = "==", double ref = 0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void indexToCoord(Rcpp::List& out, const Rcpp::IntegerVector V, int N, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void coordToIndex(Rcpp::IntegerVector& out, const Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void linRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix & rho, const Rcpp::NumericMatrix& lr, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void lltRcpp(Rcpp::NumericVector& out,const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void urtRcpp(Rcpp::NumericVector& out,const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void labRcpp(Rcpp::List & out, int nfeats, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void half2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void vector2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X, const Rcpp::IntegerVector& i, const Rcpp::IntegerVector& j, int nfeats, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void ratiosRcpp(Rcpp::NumericMatrix & out, const Rcpp::NumericMatrix & X, propr_context context=DEFAULT_GLOBAL_CONTEXT);
            
            __forceinline__
            void results2matRcpp(Rcpp::NumericMatrix & out, const Rcpp::DataFrame& results, int n, double diagonal = 0.0, propr_context context=DEFAULT_GLOBAL_CONTEXT);
        }
    }
}