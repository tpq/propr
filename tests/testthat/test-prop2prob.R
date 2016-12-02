library(propr)
context("prop2prob")

data(caneToad.counts)
rho <- perb(caneToad.counts, ivar = 15, select = 1:20)
ct1 <- caneToad.counts[1:10, 1:20]
ct2 <- caneToad.counts[11:20, 1:20]
rho1 <- perb(ct1, ivar = 15)
rho2 <- perb(ct2, ivar = 15)

dt <- prop2prob(rho1, rho2)
i <- dt$Partner[1]
j <- dt$Pair[1]

z1 <- atanh(rho1@matrix[i, j])
z2 <- atanh(rho2@matrix[i, j])

test_that("prop2prob correctly calculates probability", {

  expect_equal(
    propr:::linRcpp(rho1@matrix, rho1@logratio[])[i, j],
    z1
  )

  expect_equal(
    z1 - propr:::linRcpp(rho2@matrix, rho2@logratio[])[i, j],
    z1 - z2
  )

  r1 <- cor(rho1@logratio[, c(i, j)])[1,2]
  r2 <- cor(rho2@logratio[, c(i, j)])[1,2]
  v1 <- (1/8) * (1-r1^2) * rho1@matrix[i, j]^2 /
    ((1 - rho1@matrix[i, j]^2) * r1^2)
  v2 <- (1/8) * (1-r2^2) * rho2@matrix[i, j]^2 /
    ((1 - rho2@matrix[i, j]^2) * r2^2)
  sd <- sqrt(v1 + v2)

  expect_equal(
    propr:::linRcpp(rho1@matrix, rho1@logratio[])[j, i],
    v1
  )

  expect_equal(
    v1 + propr:::linRcpp(rho2@matrix, rho2@logratio[])[j, i],
    v1 + v2
  )

  q <- qnorm(dt$Probability[1], z1 - z2, sd)

  expect_equal(
    pnorm(q, z1 - z2, sd),
    dt$Probability[1]
  )

  expect_equal(
    pnorm(q, z1 - z2, sd),
    pnorm((q-(z1-z2))/sd)
  )

  p <- pnorm(abs(z1-z2)/sd, lower.tail = FALSE) * 2

  expect_equal(
    dt$Probability[1],
    p
  )
})

test_that("abstract correctly combines propr objects (alr)", {

  expect_equal(
    abstract(rho1, rho2)@counts,
    rho@counts
  )

  expect_equal(
    abstract(rho1, rho2)@logratio,
    rho@logratio
  )

  expect_equal(
    abstract(rho1, rho2)@matrix[i, j],
    tanh(z1 - z2)
  )
})

test_that("coordToIndex works as inverse to indexToCoord", {

  v <- c(1, 5, 10, 15, 16)
  N <- nrow(rho@matrix)
  i <- propr:::indexToCoord(v, N)

  expect_equal(
    propr:::coordToIndex(i[[1]], i[[2]], N),
    v
  )
})

test_that("abstract calls simplify correctly", {

  a1 <- abstract(rho1, rho2, dt, colBy = "Probability", cutoff = .05)
  a2 <- abstract(rho1, rho2)
  small <- dt[dt$Probability < .05,]

  expect_equal(
    a2@matrix[propr:::coordToIndex(small$Partner, small$Pair, nrow(rho1@matrix))],
    a1@matrix[a1@pairs]
  )
})
