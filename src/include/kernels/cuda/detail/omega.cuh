#pragma once

#include <Rcpp.h>

namespace propr {
    namespace detail{
        namespace cuda {
            /**
             * @brief  CUDA implementation of the “basic” Omega:  sum_{r} W(r,i)*W(r,j)
             *
             * @param  h_W    [in]  Pointer to host array of length nRows*nFeats (column-major).
             * @param  h_res  [out] Pointer to host array of length nFeats*(nFeats-1)/2.
             *                  Will be filled with the strictly lower-triangle of W^T·W.
             * @param  nRows  [in]  Number of rows in W.
             * @param  nFeats [in]  Number of columns in W.
             */
            extern "C" void
            dof_group( const double* h_W, double* h_res, int nRows, int nFeats );

            /**
             * @brief  CUDA implementation of the “weighted” Omega:
             *         For each pair (i>j), compute:
             *           Wij = W(:,i) * W(:,j)      (element-wise product)
             *           n   = sum(Wij)
             *           out = n – sum(Wij^2)/n
             *
             * @param  h_W    [in]  Pointer to host array of length nRows*nFeats (column-major).
             * @param  h_res  [out] Pointer to host array of length nFeats*(nFeats-1)/2.
             *                  Will be filled with the modified values for each (i>j).
             * @param  nRows  [in]  Number of rows in W.
             * @param  nFeats [in]  Number of columns in W.
             */
            extern "C" void
            dof_population( const double* h_W, double* h_res, int nRows, int nFeats );
        }
    }
}