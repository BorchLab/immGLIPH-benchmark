#!/usr/bin/env Rscript
# Run immGLIPH on the prepped Glanville and Huang inputs and persist a flat
# cluster-assignment TSV plus the full result list as RDS.
#
# Output files:
#   results/raw/immgliph_<dataset>_assignments.tsv  (CDR3b, cluster_id, tool, dataset)
#   results/raw/immgliph_<dataset>_full.rds         (full runGLIPH result list)
#
# Method/parameter choices match the original paper for each dataset:
#   - Glanville: method = "gliph1" (RRS local + Hamming-cutoff global +
#     connected-component clustering), reference repertoire human_v1.0_CD48
#   - Huang:     method = "gliph2" (Fisher local + Fisher global + isolated
#     clustering), reference repertoire human_v2.0_CD48
#
# Usage:
#   Rscript scripts/03_run_immgliph.R --dataset glanville
#   Rscript scripts/03_run_immgliph.R --dataset huang

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(fs)
  library(BiocParallel)
  library(immGLIPH)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = NULL,
              help = "One of {glanville, huang}."),
  make_option("--cores", type = "integer", default = 4,
              help = "BiocParallel workers. Default 4."),
  make_option("--seed", type = "integer", default = 42,
              help = "RNG seed. Default 42.")
)))

stopifnot(opt$dataset %in% c("glanville", "huang"))

set.seed(opt$seed)
if (opt$cores <= 1) {
  register(SerialParam(RNGseed = opt$seed))
} else {
  register(MulticoreParam(workers = opt$cores, RNGseed = opt$seed))
}

# On error: dump traceback so failures inside immGLIPH are diagnosable.
options(error = function() {
  cat("\n--- TRACEBACK ---\n")
  traceback(3)
  if (!interactive()) quit(status = 1)
})

cfg <- switch(opt$dataset,
  glanville = list(
    # Use gliph2 (Fisher local + isolated clustering) as the primary Glanville
    # comparison even though the published reference was produced by GLIPH1.
    # Reason: raw gliph1 with default thresholds collapses this single-cell
    # dataset into a single connected component (the published 43 groups
    # came from Glanville's post-hoc filtering on donor diversity, V-gene
    # enrichment, etc., which immGLIPH does not apply automatically). gliph2
    # uses isolated-cluster construction and produces well-separated groups
    # whose topology can be compared against the published reference. The
    # gliph1 raw run is preserved at results/raw/immgliph_glanville_gliph1_full.rds
    # for transparency and is discussed in report.qmd.
    input  = "data/glanville2017/cdr3_input.tsv",
    method = "gliph2",
    refdb  = "human_v1.0_CD48"
  ),
  huang = list(
    input  = "data/huang2020/cdr3_input.tsv",
    method = "gliph2",
    refdb  = "human_v2.0_CD48"
  )
)

input <- read_tsv(cfg$input, show_col_types = FALSE) %>%
  as.data.frame()

# Workaround for an immGLIPH 0.99.3 bug in clusterScoring.R:430-432, where
# `sample(cdr3_sequences$counts, num_members, replace = FALSE)` errors when
# num_members exceeds length(counts). Dropping the counts column sets
# counts_info = FALSE and skips the buggy clonal-expansion-score block;
# cluster membership (the only output we use for the concordance metrics)
# is unaffected. Tracked at: NEEDS-IMMGLIPH-ISSUE.
input$counts <- NULL

message(sprintf("Running immGLIPH on %s (n=%d, method=%s, refdb=%s)...",
                opt$dataset, nrow(input), cfg$method, cfg$refdb))
t0 <- Sys.time()

result <- runGLIPH(
  cdr3_sequences = input,
  method         = cfg$method,
  refdb_beta     = cfg$refdb,
  result_folder  = "",
  verbose        = TRUE
)

message(sprintf("immGLIPH finished in %.1f sec.",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Flatten cluster_list (named list of data.frames, one per cluster) into a
# tidy assignment table.
cluster_list <- result$cluster_list
if (length(cluster_list) == 0) {
  stop("immGLIPH returned no clusters for ", opt$dataset)
}

assignments <- bind_rows(lapply(names(cluster_list), function(cl) {
  tibble(
    CDR3b      = unique(cluster_list[[cl]]$CDR3b),
    cluster_id = cl
  )
})) %>%
  mutate(tool = "immgliph", dataset = opt$dataset) %>%
  distinct()

out_assign <- file.path("results/raw",
                        sprintf("immgliph_%s_assignments.tsv", opt$dataset))
out_rds    <- file.path("results/raw",
                        sprintf("immgliph_%s_full.rds", opt$dataset))
dir_create(path_dir(out_assign))
write_tsv(assignments, out_assign)
saveRDS(result, out_rds)

message(sprintf(
  "%s: %d clusters, %d (CDR3, cluster) pairs, %d unique CDR3s clustered.",
  opt$dataset,
  n_distinct(assignments$cluster_id),
  nrow(assignments),
  n_distinct(assignments$CDR3b)
))
