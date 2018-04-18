randomNum <- sample(1:1000, size = 25 * 10, replace = TRUE)
X <- matrix(randomNum, nrow = 25, ncol = 10)

identical(proprPhit(X),
          propr:::proprPhit(t(X)))

identical(proprPerb(X),
          propr:::proprPerb(t(X)))

identical(proprVLR(X),
          propr:::proprVLR(X))

identical(proprCLR(X),
          propr:::proprCLR(X))

identical(proprALR(X, 1),
          propr:::proprALR(X, 1))

set.seed(1)
a <- phit(X)
set.seed(1)
b <- propr:::phit(t(X))
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- phit(X, iter = 5, iterHow = 1)
set.seed(1)
b <- propr:::phit(t(X), iter = 5, iterHow = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- phit(X, iter = 5, iterHow = 2)
set.seed(1)
b <- propr:::phit(t(X), iter = 5, iterHow = 2)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- phit(X, iter = 5, iterHow = 1, iterSize = 5)
set.seed(1)
b <- propr:::phit(t(X), iter = 5, iterHow = 1, iterSize = 5)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- phit(X, iter = 5, iterHow = 2, iterSize = 5)
set.seed(1)
b <- propr:::phit(t(X), iter = 5, iterHow = 2, iterSize = 5)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X)
set.seed(1)
b <- propr:::perb(t(X))
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 1)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 2)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 2)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 1, iterSize = 5)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 1, iterSize = 5)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 2, iterSize = 5)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 2, iterSize = 5)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, ivar = 1)
set.seed(1)
b <- propr:::perb(t(X), ivar = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 1, ivar = 1)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 1, ivar = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 2, ivar = 1)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 2, ivar = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 1, iterSize = 3, ivar = 1)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 1, iterSize = 3, ivar = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])

set.seed(1)
a <- perb(X, iter = 5, iterHow = 2, iterSize = 3, ivar = 1)
set.seed(1)
b <- propr:::perb(t(X), iter = 5, iterHow = 2, iterSize = 3, ivar = 1)
all(a@counts - t(b@counts) == 0)
all(a@logratio - t(b@logratio) == 0)
all(a@matrix - b@matrix == 0)
identical(a@pairs[, 3:ncol(a@pairs)],
          b@pairs[, 3:ncol(a@pairs)])
