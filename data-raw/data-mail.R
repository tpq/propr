# Geelong, Torquay, Angelsea, Lorne
zips <- c("3220", "3228", "3230", "3232")

# Generate random matrix
set.seed(1)
N <- 4
M <- seq(from = 5, to = 15, length.out = N)
T <- M * rnorm(N, mean = 1, sd = 0.1)
W <- rnorm(N, mean = 10)
R <- rnorm(N, mean = 10)
F <- rep(9, N)
X <- data.frame(M, T, W, R, F)

# Label matrix with zip codes
rownames(X) <- zips
X[c("3220", "3228"), ] <- X[c("3220", "3228"), ] * 2
X <- round(X * 10)

# Make divisible by weights
X$T <- X$T - X$T %% 4
X$W <- X$W - X$W %% 5
X$R <- X$R - X$R %% 2

mail <- t(X)

devtools::use_data(mail)
