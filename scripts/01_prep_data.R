#!/usr/bin/env Rscript
# Normalize Glanville 2017 and Huang 2020 supplementary data to the common
# input schema. Source sheets were chosen to maximize overlap with each
# paper's published reference cluster output (see SOURCE.md per dataset).
#
# Output schema (`cdr3_input.tsv`, tab-separated):
#   CDR3b   TRBV   TRBJ   patient   HLA   counts
#
# Output schema (`antigen_labels.tsv`):
#   CDR3b   antigen   peptide   mhc_allele
#
# Filters applied (matching paper Methods):
#   - CDR3b starts with C, ends with F
#   - 8 <= nchar(CDR3b) <= 30
#   - drop rows missing CDR3b
#
# Usage:
#   Rscript scripts/01_prep_data.R --dataset glanville
#   Rscript scripts/01_prep_data.R --dataset huang

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(readr)
  library(fs)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = NULL,
              help = "One of {glanville, huang}.")
)))

stopifnot(opt$dataset %in% c("glanville", "huang"))

valid_cdr3 <- function(x) {
  !is.na(x) &
    str_detect(x, "^C") &
    str_detect(x, "F$") &
    nchar(x) >= 8 &
    nchar(x) <= 30
}

prep_glanville <- function() {
  # MOESM4 sheet "Sheet1" — single-cell paired-chain TCR table; this is the
  # input set the published GLIPH groups (MOESM6) were derived from. 100% of
  # the 172 reference CDR3s are present here.
  src <- "data/glanville2017/raw/41586_2017_BFnature22976_MOESM4_ESM.xlsx"
  stopifnot(file.exists(src))

  raw <- suppressMessages(suppressWarnings(
    read_excel(src, sheet = "Sheet1", .name_repair = "unique")
  ))

  raw <- raw %>%
    rename(CDR3b = CDR3beta) %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(
      CDR3b,
      TRBV    = Vbeta,
      TRBJ    = Jbeta,
      patient = as.character(Donor),
      HLA     = NA_character_,
      counts  = suppressWarnings(as.integer(BetaReads))
    ) %>%
    mutate(counts = ifelse(is.na(counts) | counts < 1, 1L, counts)) %>%
    distinct()

  # Stim is a stimulation/sort condition (MtbLys, MegaIL2, PepLib, MegaIFNg,
  # PMA), not a true antigen label, but the only per-CDR3 categorical signal
  # in MOESM4. Carried as antigen for completeness; downstream antigen
  # metrics treat it as a coarse biological category.
  labels_tbl <- raw %>%
    filter(!is.na(Stim)) %>%
    transmute(
      CDR3b,
      antigen    = Stim,
      peptide    = NA_character_,
      mhc_allele = NA_character_
    ) %>%
    distinct()

  out_dir <- "data/glanville2017"
  write_tsv(input_tbl, file.path(out_dir, "cdr3_input.tsv"))
  write_tsv(labels_tbl, file.path(out_dir, "antigen_labels.tsv"))

  message(sprintf(
    "Glanville: %d input rows, %d distinct CDR3b, %d labeled rows.",
    nrow(input_tbl), n_distinct(input_tbl$CDR3b), nrow(labels_tbl)
  ))
}

prep_huang <- function() {
  raw_dir <- "data/huang2020/raw"
  src <- file.path(raw_dir, "41587_2020_505_MOESM3_ESM.xlsx")
  stopifnot(file.exists(src))

  raw <- read_excel(src, sheet = "bulk TCR", .name_repair = "minimal") %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(
      CDR3b,
      TRBV    = Vb,
      TRBJ    = Jb,
      patient = as.character(Individual),
      HLA     = NA_character_,
      counts  = suppressWarnings(as.integer(Counts))
    ) %>%
    mutate(counts = ifelse(is.na(counts) | counts < 1, 1L, counts)) %>%
    distinct()

  known_src <- file.path(raw_dir, "41587_2020_505_MOESM4_ESM.xlsx")
  known <- read_excel(known_src, sheet = "known_CDR3", .name_repair = "minimal",
                      col_names = c("CDR3b", "TRBV", "TRBJ", "peptide"))

  labels_tbl <- input_tbl %>%
    transmute(CDR3b, antigen = "Mtb",
              peptide = NA_character_, mhc_allele = NA_character_) %>%
    rows_update(
      known %>% transmute(CDR3b, antigen = "Mtb (known epitope)",
                          peptide = peptide, mhc_allele = NA_character_),
      by = "CDR3b", unmatched = "ignore"
    ) %>%
    distinct()

  out_dir <- "data/huang2020"
  write_tsv(input_tbl, file.path(out_dir, "cdr3_input.tsv"))
  write_tsv(labels_tbl, file.path(out_dir, "antigen_labels.tsv"))

  message(sprintf(
    "Huang: %d input rows, %d distinct CDR3b, %d labeled rows.",
    nrow(input_tbl), n_distinct(input_tbl$CDR3b), nrow(labels_tbl)
  ))
}

switch(opt$dataset,
  glanville = prep_glanville(),
  huang     = prep_huang()
)
