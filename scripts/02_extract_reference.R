#!/usr/bin/env Rscript
# Extract the published reference GLIPH/GLIPH2 cluster assignments from the
# original papers' supplementary data. These are the canonical "what GLIPH
# said" outputs used as the reference vector against immGLIPH.
#
# Output schema (TSV, header row):
#   CDR3b   cluster_id   tool   dataset
#
# Output files:
#   results/raw/reference_gliph_glanville.tsv
#   results/raw/reference_gliph2_huang.tsv

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(fs)
})

extract_glanville <- function() {
  src <- "data/glanville2017/raw/41586_2017_BFnature22976_MOESM6_ESM.xlsx"
  stopifnot(file.exists(src))

  # 51 GLIPH groups, each with a space-delimited list of member CDR3s.
  scoring <- read_excel(src, sheet = "all GLIPH Group Scoring",
                        .name_repair = "minimal")

  assignments <- scoring %>%
    select(Group, members = `CDR3bs, all clones`) %>%
    filter(!is.na(Group), !is.na(members)) %>%
    mutate(CDR3b = str_split(str_squish(members), "\\s+")) %>%
    select(Group, CDR3b) %>%
    unnest(CDR3b) %>%
    filter(nzchar(CDR3b)) %>%
    transmute(
      CDR3b,
      cluster_id = paste0("gliph_", Group),
      tool = "gliph",
      dataset = "glanville"
    ) %>%
    distinct()

  out <- "results/raw/reference_gliph_glanville.tsv"
  dir_create(path_dir(out))
  write_tsv(assignments, out)

  message(sprintf(
    "Glanville GLIPH reference: %d (CDR3, cluster) pairs across %d clusters and %d unique CDR3s.",
    nrow(assignments),
    n_distinct(assignments$cluster_id),
    n_distinct(assignments$CDR3b)
  ))
}

extract_huang <- function() {
  src <- "data/huang2020/raw/41587_2020_505_MOESM5_ESM.xlsx"
  stopifnot(file.exists(src))

  members <- read_excel(src, sheet = "GLIPH_group_member",
                        .name_repair = "minimal")

  assignments_raw <- members %>%
    select(index, CDR3b = TcRb) %>%
    filter(!is.na(index), !is.na(CDR3b), nzchar(CDR3b)) %>%
    transmute(
      CDR3b,
      cluster_id = paste0("gliph2_", index),
      tool = "gliph2",
      dataset = "huang"
    ) %>%
    distinct()

  out <- "results/raw/reference_gliph2_huang.tsv"
  dir_create(path_dir(out))
  write_tsv(assignments_raw, out)

  message(sprintf(
    "Huang GLIPH2 reference (raw, all 4185 clusters): %d pairs, %d clusters, %d unique CDR3s.",
    nrow(assignments_raw),
    n_distinct(assignments_raw$cluster_id),
    n_distinct(assignments_raw$CDR3b)
  ))

  # Also extract the curated 354-cluster set (the published canonical
  # reference), filtered by the paper's criteria of: >=3 unique TCRs,
  # >=3 individuals, V-gene p < 0.05.
  filtered <- read_excel(src, sheet = "Filtered_list(354)",
                         .name_repair = "minimal")
  keep_indices <- unique(filtered$index)
  keep_ids <- paste0("gliph2_", keep_indices)

  assignments_filt <- assignments_raw %>%
    filter(cluster_id %in% keep_ids)

  out2 <- "results/raw/reference_gliph2_huang_filtered.tsv"
  write_tsv(assignments_filt, out2)

  message(sprintf(
    "Huang GLIPH2 reference (filtered 354 published groups): %d pairs, %d clusters, %d unique CDR3s.",
    nrow(assignments_filt),
    n_distinct(assignments_filt$cluster_id),
    n_distinct(assignments_filt$CDR3b)
  ))
}

extract_glanville()
extract_huang()
