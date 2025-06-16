library(testthat)
library(propr)

mat <- matrix(sample(0:10, replace = TRUE, size = 150), 10, 15)
cts <- apply(mat, 2, function(x) sum(x == 0))

df     <- data.frame(propr:::labRcpp(ncol(mat)), "Z" = propr:::ctzRcpp(mat))
df_gpu <- data.frame(propr:::labRcpp(ncol(mat)), "Z" = propr:::ctzRcpp(mat, use_gpu=TRUE))

test_that("ctzRcpp correctly counts joint zero frequency", {

  for(i in 1:nrow(df)){
    expect_equal(
      df$Z[i],
      cts[df$Partner[i]] + cts[df$Pair[i]]
    )
  }
})


test_that("ctzRcpp gpu correctly counts joint zero frequency", {
  for(i in 1:nrow(df_gpu)){
    expect_equal(
      df_gpu$Z[i],
      cts[df_gpu$Partner[i]] + cts[df_gpu$Pair[i]]
    )
  }
})
