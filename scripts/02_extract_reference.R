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

  assignments <- members %>%
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
  write_tsv(assignments, out)

  message(sprintf(
    "Huang GLIPH2 reference: %d (CDR3, cluster) pairs across %d clusters and %d unique CDR3s.",
    nrow(assignments),
    n_distinct(assignments$cluster_id),
    n_distinct(assignments$CDR3b)
  ))
}

extract_glanville()
extract_huang()
