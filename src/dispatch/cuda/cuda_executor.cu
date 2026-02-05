#include <propr/runtime/cuda_executor.hpp>

#include <propr/kernels/cuda/dispatch/backend.cuh>
#include <propr/kernels/cuda/dispatch/comparison.cuh>
#include <propr/kernels/cuda/dispatch/ctzRcpp.cuh>
#include <propr/kernels/cuda/dispatch/graflex.cuh>
#include <propr/kernels/cuda/dispatch/lr2propr.cuh>
#include <propr/kernels/cuda/dispatch/lrm.cuh>
#include <propr/kernels/cuda/dispatch/lrv.cuh>
#include <propr/kernels/cuda/dispatch/omega.cuh>

namespace propr::runtime::cuda_executor {


void wtmRcpp(double& out, const Rcpp::NumericVector& x, const Rcpp::NumericVector& w) {
    dispatch::cuda::wtmRcpp(out, x, w);
}

void wtvRcpp(double& out, const Rcpp::NumericVector& x, const Rcpp::NumericVector& w) {
    dispatch::cuda::wtvRcpp(out, x, w);
}

void corRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::corRcpp(out, X);
}

void covRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, int norm_type) {
    dispatch::cuda::covRcpp(out, X, norm_type);
}

void vlrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::vlrRcpp(out, X);
}

void clrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::clrRcpp(out, X);
}

void alrRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, int ivar) {
    dispatch::cuda::alrRcpp(out, X, ivar);
}

void symRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::symRcpp(out, X);
}

void phiRcpp(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& X, bool sym) {
    dispatch::cuda::phiRcpp(out, X, sym);
}

void rhoRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X, const Rcpp::NumericMatrix& lr, int ivar) {
    dispatch::cuda::rhoRcpp(out, X, lr, ivar);
}

void indexPairs(std::vector<int>& out, const Rcpp::NumericMatrix& X, Rcpp::String op, double ref) {
    dispatch::cuda::indexPairs(out, X, op, ref);
}

void indexToCoord(Rcpp::List& out, const Rcpp::IntegerVector V, int N) {
    dispatch::cuda::indexToCoord(out, V, N);
}

void coordToIndex(Rcpp::IntegerVector& out, const Rcpp::IntegerVector row, Rcpp::IntegerVector col, int N) {
    dispatch::cuda::coordToIndex(out, row, col, N);
}

void linRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& rho, const Rcpp::NumericMatrix& lr) {
    dispatch::cuda::linRcpp(out, rho, lr);
}

void lltRcpp(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::lltRcpp(out, X);
}

void urtRcpp(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::urtRcpp(out, X);
}

void labRcpp(Rcpp::List& out, int nfeats) {
    dispatch::cuda::labRcpp(out, nfeats);
}

void half2mat(Rcpp::NumericMatrix& out, const Rcpp::NumericVector& X) {
    dispatch::cuda::half2mat(out, X);
}

void vector2mat(
    Rcpp::NumericMatrix& out,
    const Rcpp::NumericVector& X,
    const Rcpp::IntegerVector& i,
    const Rcpp::IntegerVector& j,
    int nfeats) {
    dispatch::cuda::vector2mat(out, X, i, j, nfeats);
}

void ratiosRcpp(Rcpp::NumericMatrix& out, const Rcpp::NumericMatrix& X) {
    dispatch::cuda::ratiosRcpp(out, X);
}

void results2matRcpp(Rcpp::NumericMatrix& out, const Rcpp::DataFrame& results, int n, double diagonal) {
    dispatch::cuda::results2matRcpp(out, results, n, diagonal);
}

void lrm_basic(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y) {
    dispatch::cuda::lrm_basic(out, Y);
}

void lrm_weighted(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, Rcpp::NumericMatrix& W) {
    dispatch::cuda::lrm_weighted(out, Y, W);
}

void lrm_alpha(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, double a, Rcpp::NumericMatrix& Yfull) {
    dispatch::cuda::lrm_alpha(out, Y, a, Yfull);
}

void lrm_alpha_weighted(
    Rcpp::NumericVector& out,
    Rcpp::NumericMatrix& Y,
    Rcpp::NumericMatrix& W,
    double a,
    Rcpp::NumericMatrix& Yfull,
    Rcpp::NumericMatrix& Wfull) {
    dispatch::cuda::lrm_alpha_weighted(out, Y, W, a, Yfull, Wfull);
}

void lrv_basic(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y) {
    dispatch::cuda::lrv_basic(out, Y);
}

void lrv_weighted(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, Rcpp::NumericMatrix& W) {
    dispatch::cuda::lrv_weighted(out, Y, W);
}

void lrv_alpha(Rcpp::NumericVector& out, Rcpp::NumericMatrix& Y, double a, Rcpp::NumericMatrix& Yfull) {
    dispatch::cuda::lrv_alpha(out, Y, a, Yfull);
}

void lrv_alpha_weighted(
    Rcpp::NumericVector& out,
    Rcpp::NumericMatrix& Y,
    Rcpp::NumericMatrix& W,
    double a,
    Rcpp::NumericMatrix& Yfull,
    Rcpp::NumericMatrix& Wfull) {
    dispatch::cuda::lrv_alpha_weighted(out, Y, W, a, Yfull, Wfull);
}

void dof_global(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W) {
    dispatch::cuda::dof_global(out, W);
}

void dof_population(Rcpp::NumericVector& out, const Rcpp::NumericMatrix& W) {
    dispatch::cuda::dof_population(out, W);
}

int count_less_than(Rcpp::NumericVector& x, double cutoff) {
    return dispatch::cuda::count_less_than(x, cutoff);
}

int count_greater_than(Rcpp::NumericVector& x, double cutoff) {
    return dispatch::cuda::count_greater_than(x, cutoff);
}

int count_less_equal_than(Rcpp::NumericVector& x, double cutoff) {
    return dispatch::cuda::count_less_equal_than(x, cutoff);
}

int count_greater_equal_than(Rcpp::NumericVector& x, double cutoff) {
    return dispatch::cuda::count_greater_equal_than(x, cutoff);
}

void ctzRcpp(Rcpp::NumericVector& out, Rcpp::NumericMatrix& X) {
    dispatch::cuda::ctzRcpp(out, X);
}

void getOR(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G) {
    dispatch::cuda::getOR(out, A, G);
}

void getORperm(
    Rcpp::NumericVector& out,
    const Rcpp::IntegerMatrix& A,
    const Rcpp::IntegerMatrix& G,
    const Rcpp::IntegerVector& perm) {
    dispatch::cuda::getORperm(out, A, G, perm);
}

void permuteOR(Rcpp::NumericMatrix& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p) {
    dispatch::cuda::permuteOR(out, A, G, p);
}

void getFDR(Rcpp::List& out, double actual, const Rcpp::NumericVector& permuted) {
    dispatch::cuda::getFDR(out, actual, permuted);
}

void getG(Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& Gk) {
    dispatch::cuda::getG(out, Gk);
}

void graflex(Rcpp::NumericVector& out, const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p) {
    dispatch::cuda::graflex(out, A, Gk, p);
}

void lr2vlr(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr) {
    dispatch::cuda::lr2vlr(out, lr);
}

void lr2phi(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr) {
    dispatch::cuda::lr2phi(out, lr);
}

void lr2rho(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr) {
    dispatch::cuda::lr2rho(out, lr);
}

void lr2phs(Rcpp::NumericMatrix& out, Rcpp::NumericMatrix& lr) {
    dispatch::cuda::lr2phs(out, lr);
}

}  // namespace propr::runtime::cuda_executor
