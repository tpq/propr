#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("nvbenchr", quietly = TRUE)) {
    stop("nvbenchr is not installed. Build it from thirdparty/nvbench/R before running benchmarks.")
  }
  library(nvbenchr)
  library(propr)
})

set.seed(101)

# Axes
row_axis <- c(128L,256L, 512L)
col_axis <- c(128L,256L, 512L)
node_axis <- c(12L, 24L)
perm_axis <- c(64L, 128L)
backend_axis <- c("cpu", "gpu")

# Cached inputs so setup is outside the timed region
matrix_cache <- new.env(parent = emptyenv())
graflex_cache <- new.env(parent = emptyenv())

safe_run_or_skip <- function(state, backend, expr) {
  tryCatch(
    {
      expr
      TRUE
    },
    error = function(e) {
      state$skip(sprintf("backend=%s failed: %s", backend, conditionMessage(e)))
      FALSE
    }
  )
}

run_profiled <- function(state, backend, profile, work) {
  tag <- "r/full_walltime_sec"

  # Best-effort profiling toggles; ignore errors so benchmarks still run.
  if (backend == "gpu") try(propr:::setCudaProfile(TRUE), silent = TRUE)
  try(propr:::setHostProfile(TRUE), silent = TRUE)

  on.exit({
    if (backend == "gpu") try(propr:::setCudaProfile(FALSE), silent = TRUE)
    try(propr:::setHostProfile(FALSE), silent = TRUE)
  }, add = TRUE)

  wall <- system.time({
    state$exec(work, batched = FALSE, sync = TRUE)
  })[["elapsed"]]

  # Collect profiles and attach as summaries
  host_prof <- tryCatch(propr:::consumeHostProfile(), error = function(e) NULL)
  if (is.data.frame(host_prof) && nrow(host_prof)) {
    agg <- aggregate(ms ~ name, host_prof, sum)
    for (i in seq_len(nrow(agg))) {
      state$add_summary(paste0("host/", agg$name[i], "_sec"), agg$ms[i] / 1e3)
    }
  }
  if (backend == "gpu") {
    cuda_prof <- tryCatch(propr:::consumeCudaProfile(), error = function(e) NULL)
    if (is.data.frame(cuda_prof) && nrow(cuda_prof)) {
      agg <- aggregate(ms ~ name, cuda_prof, sum)
      for (i in seq_len(nrow(agg))) {
        state$add_summary(paste0("cuda/", agg$name[i], "_sec"), agg$ms[i] / 1e3)
      }
    }
  }

  state$add_summary(tag, wall)
}

get_matrix_inputs <- function(n_rows, n_cols) {
  key <- paste(n_rows, n_cols, sep = "x")
  if (!exists(key, envir = matrix_cache, inherits = FALSE)) {
    X <- matrix(rgamma(n_rows * n_cols, shape = 4, rate = 0.35),
                nrow = n_rows, ncol = n_cols)
    colnames(X) <- paste0("f", seq_len(n_cols))

    counts <- matrix(rpois(n_rows * n_cols, lambda = 30) + 1L,
                     nrow = n_rows, ncol = n_cols)
    zero_mat <- matrix(sample(0:10, size = n_rows * n_cols, replace = TRUE),
                       nrow = n_rows, ncol = n_cols)
    weights <- runif(n_rows, 0.25, 2)

    lr <- propr:::clrRcpp(X, use_gpu = FALSE)
    rho <- propr:::rhoRcpp(X, lr, ivar = n_cols, use_gpu = FALSE)

    llt_len <- n_cols * (n_cols - 1) / 2
    half_vec <- rnorm(llt_len)

    i_idx <- sample.int(n_cols, size = min(n_cols, max(4L, n_cols - 1L)), replace = TRUE)
    j_idx <- sample.int(n_cols, size = length(i_idx), replace = TRUE)
    vector_vals <- rnorm(length(i_idx))

    res_rows <- min(llt_len, max(8L, n_cols))
    results_df <- data.frame(
      V1 = sample.int(n_cols, size = res_rows, replace = TRUE),
      V2 = sample.int(n_cols, size = res_rows, replace = TRUE),
      V3 = runif(res_rows)
    )

    index_vec <- as.integer(seq_len(min(llt_len, 64L)))

    assign(key, list(
      X = X,
      counts = counts,
      zero_mat = zero_mat,
      weights = weights,
      lr = lr,
      rho = rho,
      half_vec = half_vec,
      i_idx = i_idx,
      j_idx = j_idx,
      vector_vals = vector_vals,
      results_df = results_df,
      index_vec = index_vec
    ), envir = matrix_cache)
  }
  get(key, envir = matrix_cache, inherits = FALSE)
}

get_graflex_inputs <- function(n_nodes, permutations) {
  key <- paste(n_nodes, permutations, sep = "x")
  if (!exists(key, envir = graflex_cache, inherits = FALSE)) {
    A <- matrix(sample(c(0L, 1L), n_nodes * n_nodes, replace = TRUE, prob = c(0.6, 0.4)),
                nrow = n_nodes, ncol = n_nodes)
    A[upper.tri(A)] <- t(A)[upper.tri(A)]
    diag(A) <- 1L
    A <- matrix(as.integer(A), nrow = n_nodes, ncol = n_nodes)

    Gk <- as.integer(sample(c(0L, 1L), n_nodes, replace = TRUE))
    G <- matrix(as.integer(tcrossprod(Gk)), nrow = n_nodes, ncol = n_nodes)

    perm <- as.integer(sample.int(n_nodes, n_nodes) - 1L) # match test indexing
    actual <- runif(1)
    permuted <- runif(permutations)

    assign(key, list(
      A = A,
      G = G,
      Gk = Gk,
      perm = perm,
      actual = actual,
      permuted = permuted
    ), envir = graflex_cache)
  }
  get(key, envir = graflex_cache, inherits = FALSE)
}

register_matrix_bench <- function(name, runner, element_count = function(n_rows, n_cols) n_rows * n_cols) {
  bench <- register_benchmark(function(state) {
    n_rows <- state$get_int64("Rows")
    n_cols <- state$get_int64("Cols")
    inputs <- get_matrix_inputs(n_rows, n_cols)
    backend <- state$get_string("Backend")

    safe_run_or_skip(state, backend, {
      state$add_element_count(element_count(n_rows, n_cols))
      run_profiled(state, backend, NULL, function(launch) {
        runner(backend, inputs, n_rows, n_cols)
      })
    })
  }, name = paste0("propr.", name))

  bench$add_int64_axis("Rows", row_axis)
  bench$add_int64_axis("Cols", col_axis)
  bench$add_string_axis("Backend", backend_axis)
  invisible(bench)
}

register_graflex_bench <- function(name, runner) {
  bench <- register_benchmark(function(state) {
    n_nodes <- state$get_int64("Nodes")
    permutations <- state$get_int64("Permutations")
    inputs <- get_graflex_inputs(n_nodes, permutations)
    backend <- state$get_string("Backend")

    safe_run_or_skip(state, backend, {
      state$add_element_count(n_nodes * (n_nodes - 1) / 2)
      run_profiled(state, backend, NULL, function(launch) {
        runner(backend, inputs, permutations)
      })
    })
  }, name = paste0("propr.", name))

  bench$add_int64_axis("Nodes", node_axis)
  bench$add_int64_axis("Permutations", perm_axis)
  bench$add_string_axis("Backend", backend_axis)
  invisible(bench)
}

register_matrix_bench("wtmRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::wtmRcpp(inputs$X[, 1], inputs$weights, use_gpu = use_gpu)
})

register_matrix_bench("wtvRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::wtvRcpp(inputs$X[, 1], inputs$weights, use_gpu = use_gpu)
})

register_matrix_bench("corRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::corRcpp(inputs$X, use_gpu = use_gpu)
})

register_matrix_bench("covRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::covRcpp(inputs$X, norm_type = 0L, use_gpu = use_gpu)
})

register_matrix_bench("vlrRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::vlrRcpp(inputs$X, use_gpu = use_gpu)
})

register_matrix_bench("clrRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::clrRcpp(inputs$X, use_gpu = use_gpu)
})

register_matrix_bench("alrRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  ivar <- n_cols
  propr:::alrRcpp(inputs$X, ivar = ivar, use_gpu = use_gpu)
})

# register_matrix_bench("symRcpp", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::symRcpp(inputs$X, use_gpu = use_gpu)
# })

register_matrix_bench("phiRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::phiRcpp(inputs$X, sym = TRUE, use_gpu = use_gpu)
})

register_matrix_bench("rhoRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::rhoRcpp(inputs$X, inputs$lr, ivar = n_cols, use_gpu = use_gpu)
})

register_matrix_bench("indexPairs", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::indexPairs(inputs$X, op = "all", ref = 0, use_gpu = use_gpu)
})

register_matrix_bench("indexToCoord", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::indexToCoord(inputs$index_vec, N = n_cols, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("coordToIndex", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::coordToIndex(inputs$i_idx, inputs$j_idx, N = n_cols, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols)

register_matrix_bench("linRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::linRcpp(inputs$rho, inputs$lr, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("lltRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::lltRcpp(inputs$rho, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("urtRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::urtRcpp(inputs$rho, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("labRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::labRcpp(n_cols, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("half2mat", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::half2mat(inputs$half_vec, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("vector2mat", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::vector2mat(inputs$vector_vals, inputs$i_idx, inputs$j_idx, nfeats = n_cols, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols)

register_matrix_bench("ratiosRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::ratiosRcpp(inputs$X, use_gpu = use_gpu)
})

# register_matrix_bench("results2matRcpp", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::results2matRcpp(inputs$results_df, n = n_cols, diagonal = 0, use_gpu = use_gpu)
# }, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("count_less_than", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  cutoff <- median(inputs$X[, 1])
  propr:::count_less_than(inputs$X[, 1], cutoff = cutoff, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_rows)

register_matrix_bench("count_greater_than", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  cutoff <- median(inputs$X[, 1])
  propr:::count_greater_than(inputs$X[, 1], cutoff = cutoff, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_rows)

register_matrix_bench("count_less_equal_than", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  cutoff <- median(inputs$X[, 1])
  propr:::count_less_equal_than(inputs$X[, 1], cutoff = cutoff, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_rows)

register_matrix_bench("count_greater_equal_than", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  cutoff <- median(inputs$X[, 1])
  propr:::count_greater_equal_than(inputs$X[, 1], cutoff = cutoff, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_rows)

register_matrix_bench("ctzRcpp", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::ctzRcpp(inputs$zero_mat, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_graflex_bench("getOR", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::getOR(inputs$A, inputs$G, use_gpu = use_gpu)
})

register_graflex_bench("getORperm", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::getORperm(inputs$A, inputs$G, inputs$perm, use_gpu = use_gpu)
})

register_graflex_bench("permuteOR", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::permuteOR(inputs$A, inputs$G, p = permutations, use_gpu = use_gpu)
})

register_graflex_bench("getFDR", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::getFDR(inputs$actual, inputs$permuted, use_gpu = use_gpu)
})

register_graflex_bench("getG", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::getG(inputs$Gk, use_gpu = use_gpu)
})

register_graflex_bench("graflex", function(backend, inputs, permutations) {
  use_gpu <- backend == "gpu"
  propr:::graflex(inputs$A, inputs$Gk, p = permutations, use_gpu = use_gpu)
})

# register_matrix_bench("lr2vlr", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::lr2vlr(inputs$lr, use_gpu = use_gpu)
# }, element_count = function(n_rows, n_cols) n_cols * n_cols)

# register_matrix_bench("lr2phi", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::lr2phi(inputs$lr, use_gpu = use_gpu)
# }, element_count = function(n_rows, n_cols) n_cols * n_cols)

# register_matrix_bench("lr2rho", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::lr2rho(inputs$lr, use_gpu = use_gpu)
# }, element_count = function(n_rows, n_cols) n_cols * n_cols)

# register_matrix_bench("lr2phs", function(backend, inputs, n_rows, n_cols) {
#   use_gpu <- backend == "gpu"
#   propr:::lr2phs(inputs$lr, use_gpu = use_gpu)
# }, element_count = function(n_rows, n_cols) n_cols * n_cols)

register_matrix_bench("lrm", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::lrm(inputs$counts, inputs$counts, weighted = FALSE, a = NA_real_,
              Yfull = inputs$counts, Wfull = inputs$counts, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("lrv", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::lrv(inputs$counts, inputs$counts, weighted = FALSE, a = NA_real_,
              Yfull = inputs$counts, Wfull = inputs$counts, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("omega", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::omega(inputs$counts, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

register_matrix_bench("Omega", function(backend, inputs, n_rows, n_cols) {
  use_gpu <- backend == "gpu"
  propr:::Omega(inputs$counts, use_gpu = use_gpu)
}, element_count = function(n_rows, n_cols) n_cols * (n_cols - 1) / 2)

run_all_benchmarks(commandArgs(trailingOnly = FALSE))
