#include <Rcpp.h>

namespace propr {
    namespace rcpp {
        namespace helpers {
            
            void printRMatrix(Rcpp::NumericMatrix mat) {
                int nrow = mat.nrow();
                int ncol = mat.ncol();
                for (int i = 0; i < nrow; ++i) {
                    for (int j = 0; j < ncol; ++j) {
                        std::cout << mat(i, j) << " ";
                    } std::cout << std::endl;
                }  std::cout << std::endl;
            }

            Rcpp::NumericMatrix pad_matrix(const Rcpp::NumericMatrix& mat,
                                           int padTop, int padBottom,
                                           int padLeft, int padRight) {
                int old_nrow = mat.nrow();
                int old_ncol = mat.ncol();
                int new_nrow = old_nrow + padTop + padBottom;
                int new_ncol = old_ncol + padLeft + padRight;
                Rcpp::NumericMatrix out(new_nrow, new_ncol);  
                for (int i = 0; i < old_nrow; ++i) {
                    for (int j = 0; j < old_ncol; ++j) {
                        out(i + padTop, j + padLeft) = mat(i, j);
                    }
                }
                return out;
            }
        }
    }
}