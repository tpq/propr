library(propr)
data(caneToad.counts)
data(caneToad.groups)

df <- caneToad.counts[,1:500]
grp <- caneToad.groups
A <- propr(df)
B <- propd(df, grp)
C <- setEmergent(B)

test_that("getNetwork matches old network functions", {

  g1a <- getNetwork(thetad.object = B, thetad.cutoff = .3)
  g1b <- plot(B, cutoff = .3)
  expect_equal(
    igraph::as_adj(g1a),
    igraph::as_adj(g1b)
  )

  g1a <- getNetwork(thetae.object = C, thetae.cutoff = .1)
  g1b <- plot(C, cutoff = .1)
  expect_equal(
    igraph::as_adj(g1a),
    igraph::as_adj(g1b)
  )

  g1a <- getNetwork(A, .85, B, .3)
  g1b <- plot(B, cutoff = .3, propr = A[">", .85])
  i <- sort(rownames(igraph::as_adj(g1a)))
  expect_equal(
    igraph::as_adj(g1a)[i,i],
    igraph::as_adj(g1b)[i,i]
  )
})
