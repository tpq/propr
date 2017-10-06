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

test_that("alpha-transformed lrv matches pseudo-lrv from alphaTheta_old", {

  expect_equal(
    propr:::lrv(as.matrix(counts), as.matrix(counts), a = .1),
    propr:::alphaTheta_old(counts, group, .1)$alrv
  )

  expect_equal(
    propd(counts, group, weighted = FALSE, alpha = .1)@theta$theta,
    propr:::alphaTheta_old(counts, group, .1)$atheta
  )

  pd <- propd(counts, group, weighted = TRUE)
  expect_equal(
    propr:::lrv(as.matrix(counts), pd@weights, weighted = TRUE, a = .1),
    propr:::alphaThetaW_old(counts, group, .1, pd@weights)$lrv
  )

  expect_equal(
    propd(counts, group, weighted = TRUE, alpha = .1)@theta$theta,
    propr:::alphaThetaW_old(counts, group, .1, pd@weights)$theta
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

test_that("lrm correctly calculates log-ratio mean", {

  expect_equal(
    propr:::calculateTheta_old(mat, grp)$lrm1,
    propr:::lrm(as.matrix(mat[1:5, ]), as.matrix(mat[1:5, ]))
  )

  expect_equal(
    propr:::calculateTheta_old(mat, grp)$lrm2,
    propr:::lrm(as.matrix(mat[6:10, ]), as.matrix(mat[6:10, ]))
  )
})

if(requireNamespace("SDMTools", quietly = TRUE)){

  x <- abs(rnorm(20))
  w <- abs(rnorm(20))

  test_that("weighted mean matches SDMTools", {

    expect_equal(
      SDMTools::wt.mean(x, w),
      propr:::wtmRcpp(x, w)
    )
  })

  test_that("weighted variance matches SDMTools", {

    expect_equal(
      SDMTools::wt.var(x, w),
      propr:::wtvRcpp(x, w)
    )
  })
}

data(mail)

test_that("lrv without weights matches vlrRcpp", {

  expect_equal(
    propr:::lltRcpp(propr:::vlrRcpp(mail[])),
    propr:::lrv(mail, mail)
  )
})

if(requireNamespace("limma", quietly = TRUE)){

  data(iris)
  keep <- iris$Species %in% c("setosa", "versicolor")
  counts <- iris[keep, 1:4] * 10
  group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

  test_that("calculateTheta matches calculateThetaW_old", {

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$theta,
      propr::calculateTheta(counts, group, weighted = TRUE)$theta
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm1,
      propr::calculateTheta(counts, group, weighted = TRUE)$lrm1
    )

    expect_equal(
      propr:::calculateThetaW_old(counts, group)$lrm2,
      propr::calculateTheta(counts, group, weighted = TRUE)$lrm2
    )

    expect_equal(
      propr::propd(counts, group, weighted = TRUE)@theta$theta,
      propr::calculateTheta(counts, group, weighted = TRUE)$theta
    )
  })
}
