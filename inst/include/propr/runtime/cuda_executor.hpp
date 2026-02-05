#pragma once

#include <Rcpp.h>
#include <vector>

namespace propr {
    namespace runtime {
        namespace cuda_executor {

            void wtmRcpp(double& out, const Rcpp::NumericVector& x, const Rcpp::NumericVector& w);
            void wtvRcpp(double& out, const Rcpp::NumericVector& x, const Rcpp::NumericVector& w);
            void corRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X);
            void covRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, int norm_type);
            void vlrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X);
            void clrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X);
            void alrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, int ivar);
            void symRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X);
            void phiRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& X, bool sym);
            void rhoRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, const Rcpp::NumericMatrix& lr, int ivar);
            void indexPairs(std::vector<int>& out, const Rcpp::NumericMatrix& X, Rcpp::String op, double ref);
            void indexToCoord(Rcpp::List& out, const Rcpp::IntegerVector V, int N);
            void coordToIndex(Rcpp::IntegerVector& out, const Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N);
            void linRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& rho, const Rcpp::NumericMatrix& lr);
            void lltRcpp(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& X);
            void urtRcpp(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& X);
            void labRcpp(Rcpp::List& out, int nfeats);
            void half2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X);
            void vector2mat(
                Rcpp::NumericMatrix& out,
                const Rcpp::NumericVector& X,
                const Rcpp::IntegerVector& i,
                const Rcpp::IntegerVector& j,
                int nfeats);
            void ratiosRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X);
            void results2matRcpp(Rcpp::NumericMatrix& out, const Rcpp::DataFrame& results, int n, double diagonal);

            void lrm_basic(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y);
            void lrm_weighted(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, Rcpp::NumericMatrix& W);
            void lrm_alpha(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, double a, Rcpp::NumericMatrix& Yfull);
            void lrm_alpha_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix& Y,
                Rcpp::NumericMatrix& W,
                double a,
                Rcpp::NumericMatrix& Yfull,
                Rcpp::NumericMatrix& Wfull);

            void lrv_basic(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y);
            void lrv_weighted(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, Rcpp::NumericMatrix& W);
            void lrv_alpha(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, double a, Rcpp::NumericMatrix& Yfull);
            void lrv_alpha_weighted(
                Rcpp::NumericVector& out,
                Rcpp::NumericMatrix& Y,
                Rcpp::NumericMatrix& W,
                double a,
                Rcpp::NumericMatrix& Yfull,
                Rcpp::NumericMatrix& Wfull);

            void dof_global(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W);
            void dof_population(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W);

            int count_less_than(Rcpp::NumericVector& x, double cutoff);
            int count_greater_than(Rcpp::NumericVector& x, double cutoff);
            int count_less_equal_than(Rcpp::NumericVector& x, double cutoff);
            int count_greater_equal_than(Rcpp::NumericVector& x, double cutoff);

            void ctzRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix& X);

            void getOR(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G);
            void getORperm(
                Rcpp::NumericVector& out,
                const Rcpp::IntegerMatrix& A,
                const Rcpp::IntegerMatrix& G,
                const Rcpp::IntegerVector& perm);
            void permuteOR(Rcpp::NumericMatrix& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p);
            void getFDR(Rcpp::List& out, double actual, const Rcpp::NumericVector& permuted);
            void getG(Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& Gk);
            void graflex(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p);

            void lr2vlr(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr);
            void lr2phi(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr);
            void lr2rho(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr);
            void lr2phs(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr);

        }  // namespace cuda_executor
    }  // namespace runtime
}  // namespace propr
