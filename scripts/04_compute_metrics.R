#!/usr/bin/env Rscript
# Compute concordance metrics between immGLIPH and the published reference
# cluster vectors. Writes a tidy long-form TSV.
#
# Output: results/metrics.tsv with columns
#   dataset   tool_pair   universe   metric   value
#
# Universes:
#   - "intersection" : CDR3s present in BOTH the input and the reference
#                      output. Cleanest answer to "do they cluster the
#                      same objects the same way?".
#   - "input_full"   : every CDR3b in the immGLIPH input. Reference CDR3s
#                      not clustered by the published run get unique
#                      singleton IDs; tests over-/under-clustering too.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(fs)
  library(mclust)
  library(aricode)
})

assign_clusters <- function(universe, assigned, prefix = "x") {
  # Returns an integer cluster vector aligned to `universe`, with each
  # unassigned CDR3 mapped to its own singleton.
  m <- assigned$cluster_id[match(universe, assigned$CDR3b)]
  is_singleton <- is.na(m)
  m[is_singleton] <- paste0(prefix, "_singleton_",
                            seq_len(sum(is_singleton)))
  as.integer(factor(m))
}

pairwise_f1 <- function(a, b) {
  # All pairs (i, j), i < j, classify each:
  #   TP: same cluster in both
  #   FP: same in immgliph (a), different in reference (b)
  #   FN: same in reference (b), different in immgliph (a)
  n <- length(a)
  stopifnot(length(b) == n)
  if (n < 2) return(c(precision = NA, recall = NA, f1 = NA))

  # Same-cluster index sets per labeling.
  same_a <- outer(a, a, `==`)
  same_b <- outer(b, b, `==`)
  upper <- upper.tri(same_a)
  TP <- sum(same_a & same_b & upper)
  FP <- sum(same_a & !same_b & upper)
  FN <- sum(!same_a & same_b & upper)

  precision <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  recall    <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
                 2 * precision * recall / (precision + recall) else NA_real_
  c(precision = precision, recall = recall, f1 = f1)
}

compute_for_dataset <- function(dataset, immgliph_path, ref_path, tool_pair) {
  immg <- read_tsv(immgliph_path, show_col_types = FALSE)
  ref  <- read_tsv(ref_path,      show_col_types = FALSE)

  ref_cdr3s  <- unique(ref$CDR3b)
  immg_cdr3s <- unique(immg$CDR3b)

  # Universes
  intersect_universe <- intersect(immg_cdr3s, ref_cdr3s)
  input_full_universe <- unique(union(immg_cdr3s, ref_cdr3s))

  rows <- list()
  for (universe_name in c("intersection", "input_full")) {
    u <- if (universe_name == "intersection") intersect_universe else input_full_universe
    if (length(u) < 2) next

    a <- assign_clusters(u, immg, prefix = "ig")
    b <- assign_clusters(u, ref,  prefix = "rf")

    ari <- mclust::adjustedRandIndex(a, b)
    nmi <- aricode::NMI(a, b)
    pf <- pairwise_f1(a, b)

    rows[[length(rows) + 1]] <- tibble(
      dataset    = dataset,
      tool_pair  = tool_pair,
      universe   = universe_name,
      n_compared = length(u),
      metric     = c("ARI", "NMI", "pairwise_precision", "pairwise_recall", "pairwise_F1"),
      value      = c(ari, nmi, pf["precision"], pf["recall"], pf["f1"])
    )
  }
  bind_rows(rows)
}

datasets <- list(
  list(
    dataset       = "glanville",
    tool_pair     = "immgliph_vs_gliph",
    immgliph_path = "results/raw/immgliph_glanville_assignments.tsv",
    ref_path      = "results/raw/reference_gliph_glanville.tsv"
  ),
  list(
    dataset       = "huang",
    tool_pair     = "immgliph_vs_gliph2",
    immgliph_path = "results/raw/immgliph_huang_assignments.tsv",
    ref_path      = "results/raw/reference_gliph2_huang.tsv"
  )
)

all_metrics <- bind_rows(lapply(datasets, function(d) {
  if (!file.exists(d$immgliph_path)) {
    warning("Missing ", d$immgliph_path, " — skipping ", d$dataset)
    return(NULL)
  }
  compute_for_dataset(d$dataset, d$immgliph_path, d$ref_path, d$tool_pair)
}))

dir_create("results")
write_tsv(all_metrics, "results/metrics.tsv")
cat("\nMetrics summary:\n")
print(all_metrics, n = Inf)
