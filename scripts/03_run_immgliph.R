#!/usr/bin/env Rscript
# Run immGLIPH on a prepped dataset and persist a flat cluster-assignment
# TSV plus the full result list as RDS.
#
# Output files (named by dataset, method, params):
#   results/raw/immgliph_<dataset>_<method>_<params>_assignments.tsv
#   results/raw/immgliph_<dataset>_<method>_<params>_full.rds
#
# `params`:
#   "default"  -- immGLIPH's runGLIPH() defaults (lcminove=c(1000,100,10),
#                  lcminp=0.01, gccutoff=NULL, sim_depth=1000)
#   "paper"    -- thresholds from the originating paper:
#                  Glanville: lcminove=c(10,10,10), lcminp=0.001, gccutoff=1
#                  Huang:     lcminove=c(10,10,10), lcminp=0.001, gccutoff=1
#                  (Huang's GLIPH2 uses Fisher rather than fold-change for
#                   motif significance; gccutoff=1 matches the global "differ
#                   at the same position by exchangeable amino acids per
#                   BLOSUM62" rule, which Hamming-1 approximates.)
#
# Usage:
#   Rscript scripts/03_run_immgliph.R --dataset glanville --method gliph2 --params default
#   Rscript scripts/03_run_immgliph.R --dataset glanville --method gliph1 --params paper
#   Rscript scripts/03_run_immgliph.R --dataset huang     --method gliph2 --params default

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(readr)
  library(fs)
  library(BiocParallel)
  library(immGLIPH)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = NULL,
              help = "One of {glanville, huang}."),
  make_option("--method", type = "character", default = NULL,
              help = "One of {gliph1, gliph2}. If NULL, uses the dataset default."),
  make_option("--params", type = "character", default = "default",
              help = "One of {default, paper}. Default 'default'."),
  make_option("--cores", type = "integer", default = 4,
              help = "BiocParallel workers."),
  make_option("--seed", type = "integer", default = 42,
              help = "RNG seed.")
)))

stopifnot(opt$dataset %in% c("glanville", "huang"),
          opt$params %in% c("default", "paper"))
if (is.null(opt$method)) {
  opt$method <- "gliph2"  # safe default; Glanville benefits from gliph1 + paper
}
stopifnot(opt$method %in% c("gliph1", "gliph2"))

set.seed(opt$seed)
if (opt$cores <= 1) {
  register(SerialParam(RNGseed = opt$seed))
} else {
  register(MulticoreParam(workers = opt$cores, RNGseed = opt$seed))
}

options(error = function() {
  cat("\n--- TRACEBACK ---\n"); traceback(3)
  if (!interactive()) quit(status = 1)
})

dataset_input <- list(
  glanville = list(input = "data/glanville2017/cdr3_input.tsv",
                   refdb = "human_v1.0_CD48"),
  huang     = list(input = "data/huang2020/cdr3_input.tsv",
                   refdb = "human_v2.0_CD48")
)[[opt$dataset]]

# Paper-matched thresholds (see header for citations).
paper_params <- list(
  lcminove      = c(10, 10, 10),
  lcminp        = 0.001,
  gccutoff      = 1,
  sim_depth     = 1000,
  kmer_mindepth = 3
)

input <- read_tsv(dataset_input$input, show_col_types = FALSE) %>%
  as.data.frame()

# Workaround for clusterScoring.R bug in 0.99.3 (fixed upstream but not yet
# in our pinned commit): drop counts to skip the buggy clonal-expansion block.
input$counts <- NULL

message(sprintf(
  "Running immGLIPH on %s (n=%d, method=%s, params=%s, refdb=%s)...",
  opt$dataset, nrow(input), opt$method, opt$params, dataset_input$refdb
))

t0 <- Sys.time()

call_args <- list(
  cdr3_sequences = input,
  method         = opt$method,
  refdb_beta     = dataset_input$refdb,
  result_folder  = "",
  verbose        = TRUE
)
if (opt$params == "paper") {
  call_args <- c(call_args, paper_params)
}

result <- do.call(runGLIPH, call_args)

message(sprintf("immGLIPH finished in %.1f sec.",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cluster_list <- result$cluster_list
if (length(cluster_list) == 0) stop("No clusters returned for ", opt$dataset)

assignments <- bind_rows(lapply(names(cluster_list), function(cl) {
  tibble(CDR3b = unique(cluster_list[[cl]]$CDR3b), cluster_id = cl)
})) %>%
  mutate(tool    = paste0("immgliph_", opt$method, "_", opt$params),
         dataset = opt$dataset) %>%
  distinct()

stem <- sprintf("immgliph_%s_%s_%s", opt$dataset, opt$method, opt$params)
out_assign <- file.path("results/raw", paste0(stem, "_assignments.tsv"))
out_rds    <- file.path("results/raw", paste0(stem, "_full.rds"))
dir_create(path_dir(out_assign))
write_tsv(assignments, out_assign)
saveRDS(result, out_rds)

message(sprintf(
  "%s: %d clusters, %d (CDR3, cluster) pairs, %d unique CDR3s clustered.",
  stem,
  n_distinct(assignments$cluster_id),
  nrow(assignments),
  n_distinct(assignments$CDR3b)
))
