#include <Rcpp.h>
#include <stdint.h>
using namespace Rcpp;

// Calculates sum(LogicalVector)
// [[Rcpp::export]]
int32_t count_if(LogicalVector x){

  int32_t counter = 0;
  for(int32_t i = 0; i < x.size(); i++){

    if(x(i) == TRUE){
      counter++;
    }
  }

  return counter;
}

// Replicates permuteTheta_old
// [[Rcpp::export]]
List pairmutate(NumericMatrix counts,
                LogicalVector group){

  // Access sample function from base R
  Environment base("package:base");
  Function sample = base["sample"];

  // Stores a vector of log-ratios for each group
  int32_t n1 = count_if(group == TRUE);
  NumericVector posij(n1);
  int32_t n2 = count_if(group != TRUE);
  NumericVector negij(n2);

  // Stores log-ratio variance for each group
  int32_t nfeats = counts.ncol();
  int32_t nsubjs = counts.nrow();
  int32_t llt = nfeats * (nfeats - 1) / 2;
  NumericVector lrv1(llt);
  NumericVector lrv2(llt);
  int32_t counter = 0;

  for(int32_t i = 1; i < nfeats; i++){
    for(int32_t j = 0; j < i; j++){

      // Randomize group labels, then retrieve log-ratios
      LogicalVector grp = sample(group);
      int32_t posct = 0;
      int32_t negct = 0;
      for(int32_t k = 0; k < nsubjs; k++){
        if(grp(k) == TRUE){
          posij(posct) = log(counts(k, i) / counts(k, j));
          posct += 1;
        }else{
          negij(negct) = log(counts(k, i) / counts(k, j));
          negct += 1;
        }
      }

      // Calculate log-ratio variance for each group
      double pvar = var(posij);
      lrv1(counter) = pvar;
      double nvar = var(negij);
      lrv2(counter) = nvar;
      counter += 1;
    }
  }

  return List::create(
    _["lrv1"] = lrv1,
    _["lrv2"] = lrv2
  );
}
