#' Convert pairwise propd results into genewise
#'
#' This function converts pairwise propd results into genewise results. The resulting
#' genewise results are direct indicators of genes being differentially expressed.
#'
#' @param propd A \code{\link{propd}} object, with FDR values from updateF. Note:
#' for the moment, only theta results with F-stats are supported. Later we will look
#' on the option to get FDR values based on permutations.
#' @param pairwise_fdr FDR threshold to consider a pairwise relationship as significant.
#' Default is 0.05.
#' @param partner_fraction Numeric in (0, 1). Fraction of genes to use as
#'   partners per batch. Default 0.017. Increase for more stable ES at the
#'   cost of compute time.
#' @param n_iter Integer. Number of independent iterations to aggregate.
#'   Default 5. Higher values give more stable results.
#' @param ... Additional arguments passed to \code{.compute_fgsea_padj_batches},
#'   such as \code{nperm}, \code{seed}, or \code{scoreType}.
#' @return A data frame with genewise results that can be used to identify differentially
#' expressed genes. It contains the following columns:
#'  - "id": gene identifier
#'  - "lfc": Log Fold Change of the gene, using the geometric mean of all genes as reference.
#'  - "lrmD": Log Ratio Mean Difference of the gene. Equivalent to the LFC, but using a subset
#'     of genes as reference (only the ones that are significantly connected to the gene).
#'  - "connectivity": number of significant pairwise relationships the gene has.
#'  - "ES": Enrichment Score of the gene.
#'  - "padj": adjusted p-value of the gene, based on the Enrichment Score.
#'
#' @rdname propdGenewise
#' @export
propdGenewise <- function(propd, pairwise_fdr = 0.05,
                          partner_fraction = 0.017,
                          n_iter = 5L, ...) {

  # for the moment it only works with theta results with F-stats,
  # but later we will look into the option to get FDR values based on permutations.
  if (!"FDR" %in% colnames(propd@results)) {
    stop("Please run updateF on the propd object to get FDR values before running propdGenewise.")
  }

  n_na_fdr <- sum(is.na(propd@results$FDR))
  if (n_na_fdr > 0) {
    stop("%d of %d pairs (%.1f%%) have NA FDR, likely due to zero counts. ",
         "Try adjusting zero-handling options (e.g., adding pseudocounts) and re-running updateF.")
  }

  # not working for alpha != NA too
  if (!is.na(propd@alpha)) {
    stop("propdGenewise currently only works for alpha = NA. Future updates may include support for other alpha values.")
  }

  # not working for more than 2 groups too
  if (length(unique(propd@group)) > 2) {
    stop("propdGenewise currently only works for 2 groups. Future updates may include support for more than 2 groups.")
  }

  # get features and number of features
  features <- colnames(propd@counts)
  nfeatures <- length(features)

  # warn for small datasets with only few genes.
  if (nfeatures < 50) {
    warning(sprintf(
      paste0("propdGenewise: only %d genes detected. \n",
             "Permutation p-values from fgsea are unreliable at this scale. \n" ,
             "Results are for testing purposes only. \n",
             "Using connnectivity is prefered. "
      ),
      nfeatures
    ))
  }

  ## ---- Build matrices needed for connectivity metrics ----
  fdr_mat <- results_to_matrix(propd@results, what = "FDR", features = features)
  theta_mat <- results_to_matrix(propd@results, what = "theta", features = features)
  adj <- (fdr_mat > 0) & (fdr_mat < pairwise_fdr)

  ## ---- Connectivity ----
  connectivity <- rowSums(adj, na.rm = TRUE)

  ## ---- GSEA-inspired enrichment score (vectorized) ----
  es_scores <- .compute_es(propd@results, features)

  ## ---- Estimate significance of ES via permutation ----
  fgsea_batches <- .compute_fgsea_padj_batches(propd, partner_fraction = 0.017,
                                               n_iter = 5L, nperm = 1000L,
                                               seed = 42L, scoreType = "pos")


  ## ---- Build lrm matrices ----
  lrm1_all <- results_to_matrix(propd@results, what = "lrm1", features = features)
  lrm2_all <- results_to_matrix(propd@results, what = "lrm2", features = features)
  # results_to_matrix returns symmetric matrices, but lrm values are directed:
  # lrm(Partner, Pair) = mean(log(x_Partner / x_Pair)), with Partner > Pair.
  # Negate the upper triangle so that mat[g, j] = mean(log(x_g / x_j)) for all g, j.
  lrm1_all[upper.tri(lrm1_all)] <- -lrm1_all[upper.tri(lrm1_all)]
  lrm2_all[upper.tri(lrm2_all)] <- -lrm2_all[upper.tri(lrm2_all)]

  ## ---- LFC (CLR-based log fold change) ----
  # lrm differences represent log fold changes; averaging across all genes as
  # reference is equivalent to using the geometric mean (CLR transformation)
  lrm_diff_full <- lrm1_all - lrm2_all
  lfc <- rowMeans(lrm_diff_full, na.rm = TRUE) / log(2)

  ## ---- lrmD (LFC using only significant partners as reference) ----
  lrm_diff_sig <- ifelse(adj, lrm_diff_full, NA) # keep only significant pairwise relationships
  lrmD <- apply(lrm_diff_sig, 1, median, na.rm = TRUE) / log(2)

  ## ---- Compile results into a data frame ----
  result <- data.frame(
    id = features,
    lfc = lfc,
    lrmD = lrmD,
    connectivity = connectivity,
    ES          = es_scores$es_pos,
    padj        = fgsea_batches$padj,
    ES_batch = fgsea_batches$ES_batch,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  ## ---- Warn if full and batch ES are poorly concordant (overall) ----
  es_cor <- cor(es_scores$es_pos, fgsea_batches$ES_batch,
                method = "spearman", use = "complete.obs")
  if (is.na(es_cor) || es_cor < 0.7) {
    warning(sprintf(
      paste0(
        "Low concordance between full ES and batch ES (Spearman rho = %.2f). ",
        "The redundancy assumption may be violated: partner subsampling in the ",
        "batch procedure may not reflect the full gene-set signal. ",
        "Treat 'padj' values with caution. Consider increasing partner_fraction ",
        "or n_iter."
      ),
      es_cor
    ))
  } else {
    message(sprintf("ES concordance (Spearman rho = %.2f): OK.", es_cor))
  }

  es_diff_threshold <- 0.5
  padj_threshold <- 0.05
  sig_mask <- result$padj < padj_threshold & !is.na(result$padj)
  relative_diff <- abs(result$ES - result$ES_batch) / (abs(result$ES) + 1e-9)
  discordant_sig <- sig_mask & (relative_diff > es_diff_threshold)

  if (any(discordant_sig)) {
    warning(sprintf(
      paste0(
        "%d significant gene(s) (padj < %.2f) show high discordance between ",
        "full ES and batch ES (relative difference > %.0f%%). ",
        "The redundancy assumption may be violated for these genes and their ",
        "'padj' values may be unreliable. ",
        "Affected genes: %s"
      ),
      sum(discordant_sig),
      padj_threshold,
      es_diff_threshold * 100,
      paste(result$id[discordant_sig], collapse = ", ")
    ))
  } else {
    message("No significant genes show high discordance between full ES and batch ES.")
  }
  return(result)
}


#' Compute GSEA-inspired enrichment scores (vectorized, integer-native)
#'
#' Pairs are ranked by theta ascending (most differentially proportional first).
#' For each gene g, its "set" is the G-1 pairs it participates in. The running
#' sum increments by 1/(G-1) on a hit and decrements by 1/(N-(G-1)) on a miss.
#'
#' Key insight: the running sum at hit j (occurring at rank r_j) has a closed form:
#'   rs(j) = j * (hit_inc + miss_inc) - r_j * miss_inc
#' The minimum always occurs just BEFORE a hit (rs(j) - hit_inc, which also
#' correctly handles the trough before the first hit) or in the tail after the
#' last hit. This gives a fully vectorized implementation via tapply.
#'
#' Partner and Pair columns in results are assumed to be 1-based integer indices
#' into features (as produced by propd), so no string matching is needed.
#'
#' @param results A data frame from propd@@results with columns Partner (int),
#'   Pair (int), theta.
#' @param features Character vector of gene names (length G).
#' @return A data frame with columns es, es_pos, es_neg (one row per gene,
#'   same order as features).
#' @keywords internal
.compute_es <- function(results, features) {

  G <- length(features)
  N <- nrow(results)

  expected_N <- G * (G - 1L) / 2L
  if (N != expected_N) {
    warning(sprintf(
      ".compute_es: expected %d pairs for %d genes but got %d. ES may be unreliable.",
      expected_N, G, N
    ))
  }

  K        <- G - 1L
  hit_inc  <- 1.0 / K
  miss_inc <- 1.0 / (N - K)

  # Rank all pairs by theta ascending (1 = most DP)
  hit_ranks <- rank(results$theta, ties.method = "first")

  # Partner and Pair are 1-based integer indices — use directly, no match() needed.
  # Build long table: each pair contributes two rows (one per gene).
  long_gene     <- c(results$Partner, results$Pair)
  long_hit_rank <- c(hit_ranks,       hit_ranks)

  # Sort by gene index, then by rank within gene
  ord           <- order(long_gene, long_hit_rank)
  long_gene     <- long_gene[ord]
  long_hit_rank <- long_hit_rank[ord]

  # Hit index j (1..K) within each gene
  j  <- as.numeric(ave(long_hit_rank, long_gene, FUN = seq_along))

  # Running sum AT hit j:  j*(hit_inc + miss_inc) - r_j * miss_inc
  rs        <- j * (hit_inc + miss_inc) - long_hit_rank * miss_inc

  # Trough just BEFORE hit j (also handles pre-first-hit descent correctly)
  rs_before <- rs - hit_inc

  # Aggregate per gene (integer keys -> tapply result indexed 1..G)
  es_pos     <- as.numeric(tapply(rs,            long_gene, max))
  last_rs    <- as.numeric(tapply(rs,            long_gene, function(x) x[length(x)]))
  last_rank  <- as.numeric(tapply(long_hit_rank, long_gene, max))
  tail_val   <- last_rs - (N - last_rank) * miss_inc
  min_before <- as.numeric(tapply(rs_before,     long_gene, min))
  es_neg     <- pmin(min_before, tail_val)

  es <- ifelse(abs(es_pos) >= abs(es_neg), es_pos, es_neg)

  data.frame(es = es, es_pos = es_pos, es_neg = es_neg, row.names = NULL)

}


#' Compute genewise ES (in batches) and permutation p-values from a propd object
#'
#' This function extends \code{\link{propdGenewise}} by adding permutation-based
#' p-values via a batched GSEA strategy. Because all gene-pair sets overlap
#' (every pair is shared by two genes), running fgsea naively would violate
#' the independence assumption of the permutation test. The solution is to
#' process genes in batches: within each batch, each focal gene is assigned a
#' disjoint random subset of partner genes, so the resulting pair sets are
#' independent by construction. fgsea is then run on this reduced, independent
#' score vector, yielding valid p-values.
#'
#' To reduce variance from the random partner assignment, the procedure is
#' repeated \code{n_iter} times with different seeds. The final ES is the mean
#' across iterations and the final p-value is the median (robust to outlier
#' iterations). Global BH correction is applied to the median p-values.
#'
#' @param propd A \code{\link{propd}} object.
#' @param partner_fraction Numeric in (0, 1). Fraction of all other genes to
#'   use as partners for each focal gene per batch. Default 0.017 (~1.7\%),
#'   giving ~300 partners at G=18080. A minimum of 100 partners is enforced.
#' @param n_iter Integer. Number of independent iterations to run and aggregate.
#'   Default 5. Higher values give more stable ES and p-values at the cost of
#'   proportionally more compute time.
#' @param nperm Integer. Number of permutations per fgseaSimple call. Default 1000.
#' @param seed Integer. Base random seed. Each iteration uses seed + i for
#'   reproducibility. Default 42.
#' @param scoreType Character. Passed to fgseaSimple. "pos" (default) tests
#'   enrichment at the low-theta end only. "std" tests both ends.
#' @return A data frame with one row per gene and columns:
#'   - "id": gene identifier
#'   - "ES_batch": median enrichment score across iterations (subsampled pairs)
#'   - "pval": median p-value across iterations
#'   - "padj": BH-adjusted median p-value (global, across all genes)
#'
#' @importFrom fgsea fgseaSimple
#' @rdname propdGenewise
#' @export
.compute_fgsea_padj_batches <- function(propd,
                                   partner_fraction = 0.017,
                                   n_iter           = 5L,
                                   nperm            = 1000L,
                                   seed             = 42L,
                                   scoreType        = "pos") {

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
  }
  if (partner_fraction <= 0 || partner_fraction >= 1) {
    stop("partner_fraction must be in (0, 1).")
  }
  n_iter <- as.integer(n_iter)
  if (n_iter < 1L) stop("n_iter must be >= 1.")

  features <- colnames(propd@counts)
  G        <- length(features)
  results  <- propd@results

  # ---- Derive max_partners ----
  max_partners <- max(100L, floor((G - 1L) * partner_fraction))
  max_partners <- as.integer(min(max_partners, G - 1L))

  batch_size <- floor(G / (1L + max_partners))
  if (batch_size < 1L) stop("Derived batch size < 1. Decrease partner_fraction.")
  n_batches <- ceiling(G / batch_size)

  message(sprintf(
    "G=%d | partner_fraction=%.3f | max_partners=%d (%.1f%%) | batch_size=%d | n_batches=%d | n_iter=%d",
    G, partner_fraction, max_partners,
    max_partners / (G - 1) * 100,
    batch_size, n_batches, n_iter
  ))

  # ---- Global score vector: -theta so high score = low theta = most DP ----
  pair_names    <- paste(results$Partner, results$Pair, sep = ",")
  global_scores <- setNames(-results$theta, pair_names)


  # ---- Run n_iter independent iterations ----
  # Each iteration returns a G x 2 matrix (ES_batch, pval) per gene
  iter_es   <- matrix(NA_real_, nrow = G, ncol = n_iter)
  iter_pval <- matrix(NA_real_, nrow = G, ncol = n_iter)

  for (iter in seq_len(n_iter)) {
    message(sprintf("Iteration %d / %d", iter, n_iter))
    iter_result <- .run_batches_once(
      results       = results,
      G             = G,
      global_scores = global_scores,
      max_partners  = max_partners,
      batch_size    = batch_size,
      n_batches     = n_batches,
      nperm         = nperm,
      seed          = seed + iter - 1L,
      scoreType     = scoreType
    )
    iter_es[,   iter] <- iter_result$ES_batch
    iter_pval[, iter] <- iter_result$pval
  }

  # ---- Aggregate across iterations ----
  # ES: median (stable summary of enrichment magnitude)
  # pval: median (robust to unlucky subsample in one iteration)
  es_median  <- apply(iter_es, 1, median, na.rm = TRUE)
  pval_median <- apply(iter_pval, 1, median, na.rm = TRUE)

  # Global BH correction on median p-values
  padj <- p.adjust(pval_median, method = "BH")

  # ---- Assemble output ----
  data.frame(
    id     = features,
    ES_batch     = es_median,
    pval   = pval_median,
    padj   = padj,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Run one iteration of the batched fgsea procedure
#'
#' Internal helper called by \code{propdGenewisePval}. Processes all genes
#' in batches, assigning disjoint partner subsets per focal gene, running
#' fgseaSimple per batch, and returning a per-gene data frame of ES and pval.
#'
#' @keywords internal
.run_batches_once <- function(results, G, global_scores, max_partners,
                              batch_size, n_batches, nperm, seed, scoreType) {

  set.seed(seed)
  all_gene_ids  <- seq_len(G)
  batch_starts  <- seq(1L, G, by = batch_size)
  batch_results <- vector("list", length(batch_starts))
  total_dropped <- 0L

  for (b in seq_along(batch_starts)) {
    batch_start <- batch_starts[b]
    batch_end   <- min(batch_start + batch_size - 1L, G)
    focal_ids   <- all_gene_ids[batch_start:batch_end]
    n_focal     <- length(focal_ids)

    available          <- sample(setdiff(all_gene_ids, focal_ids))
    effective_partners <- min(max_partners, floor(length(available) / n_focal))

    if (effective_partners < 1L) {
      warning(sprintf("Batch %d: not enough available partners, skipping.", b))
      next
    }

    pathways        <- vector("list", n_focal)
    names(pathways) <- paste0("gene_", focal_ids)

    for (i in seq_len(n_focal)) {
      partner_ids   <- available[((i - 1L) * effective_partners + 1L):
                                   (i * effective_partners)]
      higher        <- pmax(focal_ids[i], partner_ids)
      lower         <- pmin(focal_ids[i], partner_ids)
      pathways[[i]] <- paste(higher, lower, sep = ",")
    }

    batch_pair_names <- unique(unlist(pathways))
    batch_scores     <- global_scores[batch_pair_names]

    if (anyNA(batch_scores)) {
      warning(sprintf("Batch %d: %d pairs not found in results, dropping.",
                      b, sum(is.na(batch_scores))))
      batch_scores <- batch_scores[!is.na(batch_scores)]
      pathways     <- lapply(pathways, function(p) p[p %in% names(batch_scores)])
    }

    batch_scores <- sort(batch_scores, decreasing = TRUE)

    fgsea_res <- fgsea::fgseaSimple(
      pathways  = pathways,
      stats     = batch_scores,
      scoreType = scoreType,
      nperm     = nperm
    )

    fgsea_res$gene_id <- as.integer(sub("^gene_", "", fgsea_res$pathway))

    # Fill genes silently dropped by fgsea with ES=0, pval=1
    returned_pathways <- paste0("gene_", fgsea_res$gene_id)
    missing_pathways  <- setdiff(names(pathways), returned_pathways)
    if (length(missing_pathways) > 0L) {
      total_dropped <- total_dropped + length(missing_pathways)
      missing_ids   <- as.integer(sub("^gene_", "", missing_pathways))
      fgsea_res <- rbind(
        data.frame(gene_id  = fgsea_res$gene_id,
                   ES_batch = fgsea_res$ES,
                   pval     = fgsea_res$pval,
                   stringsAsFactors = FALSE),
        data.frame(gene_id  = missing_ids,
                   ES_batch = 0,
                   pval     = 1,
                   stringsAsFactors = FALSE)
      )
    } else {
      fgsea_res <- data.frame(
        gene_id  = fgsea_res$gene_id,
        ES_batch = fgsea_res$ES,
        pval     = fgsea_res$pval,
        stringsAsFactors = FALSE
      )
    }

    batch_results[[b]] <- fgsea_res
  }

  if (total_dropped > 0L) {
    warning(sprintf(
      "%d / %d genes (%.1f%%) dropped by fgseaSimple and filled with ES=0, pval=1. Consider increasing partner_fraction or nperm.",
      total_dropped, G, total_dropped / G * 100
    ))
  }

  # Combine and sort by gene_id so rows align with features order
  combined <- do.call(rbind, batch_results)
  combined[order(combined$gene_id), ]
}

