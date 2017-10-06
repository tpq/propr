#include <Rcpp.h>

using namespace Rcpp;

// Calculates sum(LogicalVector)
// [[Rcpp::export]]
int count_if(LogicalVector x){

  int counter = 0;
  for(int i = 0; i < x.size(); i++){

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
  int n1 = count_if(group == TRUE);
  NumericVector posij(n1);
  int n2 = count_if(group != TRUE);
  NumericVector negij(n2);

  // Stores log-ratio variance for each group
  int nfeats = counts.ncol();
  int nsubjs = counts.nrow();
  int llt = nfeats * (nfeats - 1) / 2;
  NumericVector lrv1(llt);
  NumericVector lrv2(llt);
  int counter = 0;

  for(int i = 1; i < nfeats; i++){
    for(int j = 0; j < i; j++){

      // Randomize group labels, then retrieve log-ratios
      LogicalVector grp = sample(group);
      int posct = 0;
      int negct = 0;
      for(int k = 0; k < nsubjs; k++){
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
