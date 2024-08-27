library(testthat)
library(propr)

# define data matrix
set.seed(123)
N <- 100
a <- seq(from = 5, to = 15, length.out = N)
b <- a * rnorm(N, mean = 1, sd = 0.1)
c <- rnorm(N, mean = 10)
d <- rnorm(N, mean = 10)
e <- rep(10, N)
X <- data.frame(a, b, c, d, e)

