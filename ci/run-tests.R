#!/usr/bin/env Rscript
options(warn = 1)

# --- Stable absolute paths (survive testthat's wd changes) ---
start_wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
lib <- Sys.getenv("R_LIBS_USER", file.path(start_wd, ".r-lib"))
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

# --- Repos: include Bioconductor so limma installs ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
options(repos = BiocManager::repositories())

# --- Ensure runner tooling ---
ensure <- function(pkgs) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) install.packages(need)
}
ensure(c("remotes","testthat","xml2"))

message("Installing package dependencies (cached in .r-lib)...")
remotes::install_deps(dependencies = TRUE, upgrade = "never")

message("Installing 'propr' from local source ...")
# <<<<< remotes::install_local(start_wd, upgrade = "never", force = TRUE)
# install.packages(start_wd, repos = NULL, type = "source")



# --- Absolute path for JUnit output ---
junit_dir  <- file.path(start_wd, "ci_artifacts")
junit_file <- file.path(junit_dir, "junit.xml")
dir.create(junit_dir, showWarnings = FALSE, recursive = TRUE)

message("Running tests (stop_on_failure=TRUE)...")
library(testthat)

reporters <- list(
  JunitReporter$new(file = junit_file),
  SummaryReporter$new()
)

test_local(
  reporter = MultiReporter$new(reporters),
  stop_on_failure = TRUE
)

message("\n✅ All tests passed.")
