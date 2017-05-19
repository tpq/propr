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
