#!/usr/bin/env Rscript
# Apply paper-style post-hoc filtering to an immGLIPH run.
#
# Glanville 2017 published 43 GLIPH groups after filtering raw clusters with:
#   * >=3 donors per cluster
#   * >=4 unique clones (CDR3s) per cluster
#   * V-gene enrichment   p < 0.05
#   * HLA enrichment      p < 0.05  (skipped when no HLA on input)
#   * CDR3-length         p < 0.05
#
# Huang 2020's curated 354 GLIPH2 set used:
#   * >=3 unique TCRs per cluster
#   * >=3 individuals per cluster
#   * V-gene enrichment   p < 0.05
#
# Inputs:
#   --full-rds  : the full result RDS produced by 03_run_immgliph.R
#   --input-tsv : the matching cdr3_input.tsv (for per-cluster donor counts)
#   --filter    : "glanville" or "huang"
#   --out       : output assignments TSV path
#
# Output: TSV with CDR3b, cluster_id, tool, dataset (matches 03_run schema).

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(readr)
  library(fs)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--full-rds",  type = "character", default = NULL),
  make_option("--input-tsv", type = "character", default = NULL),
  make_option("--filter",    type = "character", default = NULL,
              help = "One of {glanville, huang}."),
  make_option("--out",       type = "character", default = NULL)
)))

stopifnot(!is.null(opt$`full-rds`), !is.null(opt$`input-tsv`),
          opt$filter %in% c("glanville", "huang"), !is.null(opt$out))

result <- readRDS(opt$`full-rds`)
input  <- read_tsv(opt$`input-tsv`, show_col_types = FALSE)

cluster_list <- result$cluster_list
cluster_props <- result$cluster_properties

if (length(cluster_list) == 0) {
  warning("No clusters in input RDS; writing empty assignments.")
  dir_create(path_dir(opt$out))
  write_tsv(tibble(CDR3b = character(), cluster_id = character(),
                   tool = character(), dataset = character()), opt$out)
  quit(status = 0)
}

# Build per-cluster summary: cluster_id, n_unique_cdr3, n_donors, scores.
cluster_summary <- tibble(
  cluster_id = names(cluster_list),
  n_unique_cdr3 = vapply(cluster_list, function(d) length(unique(d$CDR3b)),
                         integer(1)),
  n_donors = vapply(cluster_list, function(d) {
    if ("patient" %in% names(d)) length(unique(d$patient[!is.na(d$patient)]))
    else NA_integer_
  }, integer(1))
) %>%
  left_join(cluster_props %>%
              tibble::as_tibble(rownames = NULL) %>%
              tibble::rowid_to_column("row_idx") %>%
              mutate(cluster_id = names(cluster_list)[row_idx]) %>%
              select(-row_idx),
            by = "cluster_id")

# Apply filter.
filt_clusters <- switch(opt$filter,
  glanville = cluster_summary %>%
    filter(n_unique_cdr3 >= 4,
           is.na(n_donors) | n_donors >= 3,
           vgene.score < 0.05 | is.na(vgene.score),
           cdr3.length.score < 0.05 | is.na(cdr3.length.score)),
  huang = cluster_summary %>%
    filter(n_unique_cdr3 >= 3,
           is.na(n_donors) | n_donors >= 3,
           vgene.score < 0.05 | is.na(vgene.score))
)

dataset <- if (grepl("glanville", opt$out, ignore.case = TRUE)) "glanville" else "huang"
tool_tag <- sub("_assignments\\.tsv$", "", basename(opt$out))

assignments <- bind_rows(lapply(filt_clusters$cluster_id, function(cl) {
  tibble(CDR3b = unique(cluster_list[[cl]]$CDR3b), cluster_id = cl)
})) %>%
  mutate(tool = tool_tag, dataset = dataset) %>%
  distinct()

dir_create(path_dir(opt$out))
write_tsv(assignments, opt$out)

message(sprintf(
  "%s: %d/%d clusters pass %s filter. Writing %d (CDR3, cluster) pairs covering %d unique CDR3s.",
  basename(opt$out),
  nrow(filt_clusters), nrow(cluster_summary), opt$filter,
  nrow(assignments), n_distinct(assignments$CDR3b)
))
