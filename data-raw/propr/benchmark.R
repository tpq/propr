# Set number of features in counts matrix
i <- 2000
X <- t(data.frame("a" = sample(1:i), "b" = sample(1:i), "c" = sample(1:i),
                  "d" = sample(1:i), "e" = sample(1:i), "f" = sample(1:i)))

# Benchmark performance for several calculations
microbenchmark::microbenchmark(propr:::proprPhit(X), proprPhit(X))
microbenchmark::microbenchmark(propr:::proprPerb(X), proprPerb(X))
microbenchmark::microbenchmark(propr:::proprVLR(X), proprVLR(X))

# Benchmark performance for Rcpp calculations
microbenchmark::microbenchmark(propr:::proprVLR(X), vlrRcpp(X))
microbenchmark::microbenchmark(propr:::proprCLR(X), clrRcpp(X))
microbenchmark::microbenchmark(propr:::proprALR(X, ivar = 5), alrRcpp(X, ivar = 5))
microbenchmark::microbenchmark(propr:::proprSym(propr:::proprPhit(X)), symRcpp(phiRcpp(X)))
microbenchmark::microbenchmark(propr:::proprPhit(X), phiRcpp(X))
microbenchmark::microbenchmark(propr:::proprPerb(X, ivar = 5), rhoRcpp(X, ivar = 5))
microbenchmark::microbenchmark(propr:::phit(X), phit(X))
