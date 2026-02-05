#include <stdexcept>

#include <propr/runtime/cuda_executor.hpp>

[[noreturn]] 
static void throw_unavailable() {
    throw std::runtime_error("CUDA backend is not compiled into this build.");
}

namespace propr::runtime::cuda_executor {

// backend ops
void wtmRcpp(double&, const Rcpp::NumericVector&, const Rcpp::NumericVector&) { throw_unavailable(); }
void wtvRcpp(double&, const Rcpp::NumericVector&, const Rcpp::NumericVector&) { throw_unavailable(); }
void corRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void covRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&, int) { throw_unavailable(); }
void vlrRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void clrRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void alrRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&, int) { throw_unavailable(); }
void symRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void phiRcpp(Rcpp::NumericMatrix&, Rcpp::NumericMatrix&, bool) { throw_unavailable(); }
void rhoRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&, int) { throw_unavailable(); }
void indexPairs(std::vector<int>&, const Rcpp::NumericMatrix&, Rcpp::String, double) { throw_unavailable(); }
void indexToCoord(Rcpp::List&, const Rcpp::IntegerVector, int) { throw_unavailable(); }
void coordToIndex(Rcpp::IntegerVector&, const Rcpp::IntegerVector, Rcpp::IntegerVector, int) { throw_unavailable(); }
void linRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void lltRcpp(Rcpp::NumericVector&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void urtRcpp(Rcpp::NumericVector&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void labRcpp(Rcpp::List&, int) { throw_unavailable(); }
void half2mat(Rcpp::NumericMatrix&, const Rcpp::NumericVector&) { throw_unavailable(); }
void vector2mat(Rcpp::NumericMatrix&, const Rcpp::NumericVector&, const Rcpp::IntegerVector&, const Rcpp::IntegerVector&, int) { throw_unavailable(); }
void ratiosRcpp(Rcpp::NumericMatrix&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void results2matRcpp(Rcpp::NumericMatrix&, const Rcpp::DataFrame&, int, double) { throw_unavailable(); }

// lrm ops
void lrm_basic(Rcpp::NumericVector&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrm_weighted(Rcpp::NumericVector&, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrm_alpha(Rcpp::NumericVector&, Rcpp::NumericMatrix&, double, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrm_alpha_weighted(Rcpp::NumericVector&, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&, double, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }

// lrv ops
void lrv_basic(Rcpp::NumericVector&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrv_weighted(Rcpp::NumericVector&, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrv_alpha(Rcpp::NumericVector&, Rcpp::NumericMatrix&, double, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lrv_alpha_weighted(Rcpp::NumericVector&, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&, double, Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }

// omega ops
void dof_global(Rcpp::NumericVector&, const Rcpp::NumericMatrix&) { throw_unavailable(); }
void dof_population(Rcpp::NumericVector&, const Rcpp::NumericMatrix&) { throw_unavailable(); }

// comparison ops
int count_less_than(Rcpp::NumericVector&, double) { throw_unavailable(); }
int count_greater_than(Rcpp::NumericVector&, double) { throw_unavailable(); }
int count_less_equal_than(Rcpp::NumericVector&, double) { throw_unavailable(); }
int count_greater_equal_than(Rcpp::NumericVector&, double) { throw_unavailable(); }

// ctzRcpp ops
void ctzRcpp(Rcpp::NumericVector&, Rcpp::NumericMatrix&) { throw_unavailable(); }

// graflex ops
void getOR(Rcpp::NumericVector&, const Rcpp::IntegerMatrix&, const Rcpp::IntegerMatrix&) { throw_unavailable(); }
void getORperm(Rcpp::NumericVector&, const Rcpp::IntegerMatrix&, const Rcpp::IntegerMatrix&, const Rcpp::IntegerVector&) { throw_unavailable(); }
void permuteOR(Rcpp::NumericMatrix&, const Rcpp::IntegerMatrix&, const Rcpp::IntegerMatrix&, int) { throw_unavailable(); }
void getFDR(Rcpp::List&, double, const Rcpp::NumericVector&) { throw_unavailable(); }
void getG(Rcpp::IntegerMatrix&, const Rcpp::IntegerVector&) { throw_unavailable(); }
void graflex(Rcpp::NumericVector&, const Rcpp::IntegerMatrix&, const Rcpp::IntegerVector&, int) { throw_unavailable(); }

// lr2propr ops
void lr2vlr(Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lr2phi(Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lr2rho(Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }
void lr2phs(Rcpp::NumericMatrix&, Rcpp::NumericMatrix&) { throw_unavailable(); }

}  // namespace propr::runtime::cuda_executor
