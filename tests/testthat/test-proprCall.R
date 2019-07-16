library(propr)
library(parallel)

describe("updateCutoffs.propr()", {
  num_rows <- 5
  num_cols <- 10

  # Make count-like data
  counts <- rnbinom(num_rows * num_cols, size = 0.1, mu = 10) + 1
  dat <- matrix(counts, nrow = num_rows)

  pr <- propr(dat, "phs", p = 20)

  oldResult <- propr:::updateCutoffs_old.propr(pr)
  singleCoreResult <- updateCutoffs.propr(pr, ncores = 1)
  multiCoreResult <- updateCutoffs.propr(pr, ncores = 2)

  describe("with one core", {
    it("matches the old version of the function", {
      expect_equal(singleCoreResult@fdr, oldResult@fdr)
    })
  })

  describe("with multiple cores", {
    it("matches the single core version of the function", {
      expect_equal(multiCoreResult@fdr, singleCoreResult@fdr)
    })

    it("matches the old version of the function", {
      expect_equal(multiCoreResult@fdr, oldResult@fdr)
    })

    describe("warning messages", {
      describe("with multiple cores selected", {
        it("does not warn user if parallel package is attached", {
          expect_silent(updateCutoffs.propr(pr, ncores = 2))
        })

        it("warns user if parallel package is not attached", {
          # Unload the package
          unloadNamespace("parallel")

          expect_message(updateCutoffs.propr(pr, ncores = 2))

          # Reattach package
          library("parallel")
        })
      })

      describe("with single core selected", {
        it("does not warn if parallel package is not attached", {
          # Unload the package
          unloadNamespace("parallel")

          expect_silent(updateCutoffs.propr(pr, ncores = 1))

          # Reattach package
          library("parallel")
        })
      })
    })
  })
})
