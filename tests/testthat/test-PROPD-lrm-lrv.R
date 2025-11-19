library(testthat)
library(propr)

read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){
  ext <- tolower(tail(strsplit(basename(file), split = "\\.")[[1]], 1))
  if (ext == "tsv" || ext == "txt") {
    separator <- "\t"
  } else if (ext == "csv") {
    separator <- ","
  } else {
    stop(paste("Unknown separator for", ext))
  }
  read.delim(
    file,
    sep         = separator,
    header      = header,
    row.names   = row.names,
    check.names = check.names
  )
}

counts <- read_delim_flexible(
  test_path("..", "..", "data-raw", "data", "all_counts.csv"),
  header      = TRUE,
  check.names = FALSE
)

samplesheet <- read_delim_flexible(
  test_path("..", "..", "data-raw", "data", "samplesheet.csv")
)

group <- as.vector(samplesheet[,"condition"])
group <- as.character(group)

if (length(group) != nrow(counts)) stop("Error when parsing group")
if (length(unique(group)) != 2)    stop("Only two groups are allowed for contrast")

test_that("GPU and CPU propd results are numerically equivalent", {
  ## --- CPU run ---
  options(propr.use_gpu = FALSE)
  pd_cpu <- propd(
    counts,
    group    = group,
    alpha    = 0.5,
    weighted = FALSE
  )

  ## --- GPU run ---
  options(propr.use_gpu = TRUE)
  pd_gpu <- propd(
    counts,
    group    = group,
    alpha    = 0.5,
    weighted = FALSE
  )

  expect_s4_class(pd_cpu, "propr")
  expect_s4_class(pd_gpu, "propr")

  expect_true(all(c("results") %in% slotNames(pd_cpu)))
  expect_true(all(c("results") %in% slotNames(pd_gpu)))

  cpu_res <- pd_cpu@results
  gpu_res <- pd_gpu@results

  expect_equal(dim(cpu_res), dim(gpu_res))
  expect_equal(colnames(cpu_res), colnames(gpu_res))

  ## ensure row ordering is identical
  expect_equal(cpu_res$Partner, gpu_res$Partner)
  expect_equal(cpu_res$Pair,    gpu_res$Pair)

  ## --- Numeric equality with tolerance ---
  # exclude these columns from numeric comparison
  non_num_cols <- c("Partner", "Pair")

  num_cols <- setdiff(colnames(cpu_res), non_num_cols)

  # compare columns with a small tolerance 
  expect_equal(
    as.matrix(cpu_res[, num_cols]),
    as.matrix(gpu_res[, num_cols]),
    tolerance = 1e-3,
    scale     = 1
  )
})
