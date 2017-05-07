library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

new <- propriety:::calculateTheta(counts, group)
old <- propriety:::calculateTheta_old(counts, group)

test_that("fast calculateTheta matches slow calculateTheta", {

  expect_equal(
    new[, "theta"],
    old[, "theta"]
  )
})

set.seed(1)
newp <- propriety:::permuteTheta(counts, group, p = 3)
set.seed(1)
oldp <- propriety:::permuteTheta_old(counts, group, p = 3)

test_that("fast permuteTheta matches slow permuteTheta", {

  expect_equal(
    as.vector(as.matrix(newp)),
    as.vector(oldp)
  )
})

new$theta <- new$theta^5
newp <- newp^6
old$theta <- old$theta^5
oldp <- oldp^6

test_that("new FDR matches old FDR", {

  expect_equal(
    as.vector(as.matrix(propriety:::calculateFDR(new, newp))),
    as.vector(propriety:::calculateFDR_old(old, oldp))
  )
})
