.nvtx_enabled <- function() {
  isTRUE(getOption("propr.enable_nvtx", FALSE)) && requireNamespace("nvtxR", quietly = TRUE)
}

NVTX_PUSH <- function(label, color = 0L) {
  if (.nvtx_enabled()) nvtxR::nvtx_push_range(label, as.integer(color))
  invisible(NULL)
}
NVTX_POP <- function() {
  if (.nvtx_enabled()) nvtxR::nvtx_pop_range()
  invisible(NULL)
}

WITH_NVTX <- function(label, color = 0L, expr) {
  if (.nvtx_enabled()) {
    nvtxR::nvtx_push_range(label, as.integer(color))
    on.exit(nvtxR::nvtx_pop_range(), add = TRUE)
  }
  eval.parent(substitute(expr))
}