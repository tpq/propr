# Set number of features in counts matrix
i <- 8000
counts <- data.frame("a" = sample(1:i), "b" = sample(1:i), "c" = sample(1:i),
                     "d" = sample(1:i), "e" = sample(1:i), "f" = sample(1:i))
counts <- t(counts)
X <- counts[, 1:2000]

# Benchmark performance for several calculations
microbenchmark::microbenchmark(proprPhit(X), propr:::proprPhit(X))
microbenchmark::microbenchmark(proprPerb(X), propr:::proprPerb(X))
microbenchmark::microbenchmark(proprVLR(X), propr:::proprVLR(X))
