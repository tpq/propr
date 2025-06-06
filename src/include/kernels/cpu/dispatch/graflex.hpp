#include <Rcpp.h>
#include <numeric>

namespace propr {
    namespace dispatch {
        namespace cpu {
            void getOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, Rcpp::NumericVector& out);
            void getORperm(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, const Rcpp::IntegerVector& perm, Rcpp::NumericVector& out);
            void permuteOR(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerMatrix& G, int p, Rcpp::NumericMatrix& out);
            void getFDR(double actual, const Rcpp::NumericVector& permuted, Rcpp::List& out);
            void getG(const Rcpp::IntegerVector& Gk, Rcpp::IntegerMatrix& out);
            void graflex(const Rcpp::IntegerMatrix& A, const Rcpp::IntegerVector& Gk, int p, Rcpp::NumericVector& out);
        }
    }
}