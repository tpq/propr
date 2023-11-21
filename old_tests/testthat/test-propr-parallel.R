library(propr)
library(parallel)

describe("updateCutoffs.propr()", {

  num_rows <- 5
  num_cols <- 10

  # Make count-like data
  counts <- rnbinom(num_rows * num_cols, size = 0.1, mu = 10) + 1
  dat <- matrix(counts, nrow = num_rows)

  pr <- propr(dat, "phs", p = 20)

  singleCoreResult <- updateCutoffs(pr, ncores = 1)
  multiCoreResult <- updateCutoffs(pr, ncores = 2)

  describe("with multiple cores", {
    it("matches the single core version of the function", {
      expect_equal(multiCoreResult@fdr, singleCoreResult@fdr)
    })
  })
})
