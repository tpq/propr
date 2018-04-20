library(propr)

data(iris)
keep <- iris$Species %in% c("setosa", "versicolor")
counts <- iris[keep, 1:4] * 10
group <- ifelse(iris[keep, "Species"] == "setosa", "A", "B")

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
