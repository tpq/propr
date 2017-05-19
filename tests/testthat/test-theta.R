library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

new <- propr:::calculateTheta(counts, group)
old <- propr:::calculateTheta_old(counts, group)

test_that("fast calculateTheta matches slow calculateTheta", {

  expect_equal(
    new[, "theta"],
    old[, "theta"]
  )
})

set.seed(1)
newp <- propr:::permuteTheta(counts, group, p = 3)
set.seed(1)
oldp <- propr:::permuteTheta_old(counts, group, p = 3)

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
    as.vector(as.matrix(propr:::calculateFDR(new, newp))),
    as.vector(propr:::calculateFDR_old(old, oldp))
  )
})

test_that("boxRcpp matches pseudo-lrv from alphaTheta_old", {

  expect_equal(
    propr:::boxRcpp(as.matrix(counts), .1),
    propr:::alphaTheta_old(counts, group, .1)$alrv
  )
})

test_that("fast calculateTheta matches slow alphaTheta", {

  expect_equal(
    as.vector(propr:::calculateTheta(counts, group, .1)$theta),
    as.vector(propr:::alphaTheta_old(counts, group, .1)$atheta)
  )
})

mat <- matrix(sample(0:10, replace = TRUE, size = 150), 10, 15)
cts <- apply(mat, 2, function(x) sum(x == 0))

df <- data.frame(propr:::labRcpp(ncol(mat)),
                 "Z" = propr:::ctzRcpp(mat))

test_that("ctzRcpp correctly counts joint zero frequency", {

  for(i in 1:nrow(df)){
    expect_equal(
      df$Z[i],
      cts[df$Partner[i]] + cts[df$Pair[i]]
    )
  }
})

mat <- as.data.frame(matrix(sample(1:10, replace = TRUE, size = 150), 10, 15))
grp <- c(rep("A", 5), rep("B", 5))

test_that("lrmRcpp correctly calculates log-ratio mean", {

  expect_equal(
    propr:::calculateTheta_old(mat, grp)$lrm1,
    propr:::lrmRcpp(as.matrix(mat[1:5, ]))
  )

  expect_equal(
    propr:::calculateTheta_old(mat, grp)$lrm2,
    propr:::lrmRcpp(as.matrix(mat[6:10, ]))
  )
})
