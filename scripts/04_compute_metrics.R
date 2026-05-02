#!/usr/bin/env Rscript
# Compute concordance metrics between every immGLIPH cluster vector and its
# matching published reference. Writes a tidy long-form TSV.
#
# Output: results/metrics.tsv with columns
#   dataset   comparison   universe   n_compared   metric   value
#
# Universes:
#   intersection : CDR3s in BOTH the immGLIPH input and the reference output
#                  (cleanest "do they cluster the same objects the same way")
#   input_full   : every CDR3b in either input or reference (singletons fill)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(fs)
  library(mclust)
  library(aricode)
})

assign_clusters <- function(universe, assigned, prefix = "x") {
  m <- assigned$cluster_id[match(universe, assigned$CDR3b)]
  is_singleton <- is.na(m)
  m[is_singleton] <- paste0(prefix, "_singleton_",
                            seq_len(sum(is_singleton)))
  as.integer(factor(m))
}

pairwise_f1 <- function(a, b) {
  n <- length(a)
  stopifnot(length(b) == n)
  if (n < 2) return(c(precision = NA, recall = NA, f1 = NA))
  same_a <- outer(a, a, `==`)
  same_b <- outer(b, b, `==`)
  upper <- upper.tri(same_a)
  TP <- sum(same_a & same_b & upper)
  FP <- sum(same_a & !same_b & upper)
  FN <- sum(!same_a & same_b & upper)
  precision <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  recall    <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
          2 * precision * recall / (precision + recall) else NA_real_
  c(precision = precision, recall = recall, f1 = f1)
}

compute_for_pair <- function(dataset, comparison, immgliph_path, ref_path) {
  if (!file.exists(immgliph_path) || !file.exists(ref_path)) {
    warning("Missing file(s) for ", comparison, " — skipping.")
    return(NULL)
  }
  immg <- read_tsv(immgliph_path, show_col_types = FALSE)
  ref  <- read_tsv(ref_path,      show_col_types = FALSE)

  intersect_universe  <- intersect(unique(immg$CDR3b), unique(ref$CDR3b))
  input_full_universe <- unique(union(unique(immg$CDR3b), unique(ref$CDR3b)))

  rows <- list()
  for (universe_name in c("intersection", "input_full")) {
    u <- if (universe_name == "intersection") intersect_universe else input_full_universe
    if (length(u) < 2) next
    a <- assign_clusters(u, immg, prefix = "ig")
    b <- assign_clusters(u, ref,  prefix = "rf")
    pf <- pairwise_f1(a, b)
    rows[[length(rows) + 1]] <- tibble(
      dataset    = dataset,
      comparison = comparison,
      universe   = universe_name,
      n_compared = length(u),
      metric     = c("ARI", "NMI", "pairwise_precision", "pairwise_recall", "pairwise_F1"),
      value      = c(mclust::adjustedRandIndex(a, b),
                     aricode::NMI(a, b),
                     pf["precision"], pf["recall"], pf["f1"])
    )
  }
  bind_rows(rows)
}

raw <- "results/raw"
ref_glanville      <- file.path(raw, "reference_gliph_glanville.tsv")
ref_huang_raw      <- file.path(raw, "reference_gliph2_huang.tsv")
ref_huang_filtered <- file.path(raw, "reference_gliph2_huang_filtered.tsv")

comparisons <- list(
  # Glanville: published reference is the curated 43 GLIPH groups.
  list(dataset = "glanville", comparison = "gliph2_default vs GLIPH",
       immg = file.path(raw, "immgliph_glanville_gliph2_default_assignments.tsv"),
       ref  = ref_glanville),
  list(dataset = "glanville", comparison = "gliph1_paper vs GLIPH",
       immg = file.path(raw, "immgliph_glanville_gliph1_paper_assignments.tsv"),
       ref  = ref_glanville),
  list(dataset = "glanville", comparison = "gliph1_paper+filtered vs GLIPH",
       immg = file.path(raw, "immgliph_glanville_gliph1_paper_filtered_assignments.tsv"),
       ref  = ref_glanville),

  # Huang: two reference sets, raw 4185 vs curated 354.
  list(dataset = "huang", comparison = "gliph2_default vs GLIPH2 (raw 4185)",
       immg = file.path(raw, "immgliph_huang_gliph2_default_assignments.tsv"),
       ref  = ref_huang_raw),
  list(dataset = "huang", comparison = "gliph2_default vs GLIPH2 (filtered 354)",
       immg = file.path(raw, "immgliph_huang_gliph2_default_assignments.tsv"),
       ref  = ref_huang_filtered),
  list(dataset = "huang", comparison = "gliph2_paper vs GLIPH2 (filtered 354)",
       immg = file.path(raw, "immgliph_huang_gliph2_paper_assignments.tsv"),
       ref  = ref_huang_filtered),
  list(dataset = "huang", comparison = "gliph2_paper+filtered vs GLIPH2 (filtered 354)",
       immg = file.path(raw, "immgliph_huang_gliph2_paper_filtered_assignments.tsv"),
       ref  = ref_huang_filtered)
)

all_metrics <- bind_rows(lapply(comparisons, function(c) {
  compute_for_pair(c$dataset, c$comparison, c$immg, c$ref)
}))

dir_create("results")
write_tsv(all_metrics, "results/metrics.tsv")
cat("\nMetrics summary:\n")
print(all_metrics, n = Inf)
