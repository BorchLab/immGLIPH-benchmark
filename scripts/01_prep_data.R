#!/usr/bin/env Rscript
# Normalize Glanville 2017 and Huang 2020 supplementary data to the common
# input schema consumed by every downstream tool.
#
# Output schema (`cdr3_input.tsv`, tab-separated, header row):
#   CDR3b   TRBV   TRBJ   patient   HLA   counts
#
# Output schema (`antigen_labels.tsv`):
#   CDR3b   antigen   peptide   mhc_allele
#
# Filters applied (matching what each paper used):
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
  raw_dir <- "data/glanville2017/raw"
  src <- file.path(raw_dir, "41586_2017_BFnature22976_MOESM2_ESM.xlsx")
  stopifnot(file.exists(src))

  # MOESM2 sheet "Raw" — the canonical TCR table from the paper.
  # Headers: DataSource, ID, subject, HLA, antigen-species, antigen,
  #          peptide, Reads, TCRBV, TCRBJ, CDR-H1, CDR-H2, CDR3b, ...
  raw <- read_excel(src, sheet = "Raw", .name_repair = "minimal")

  # Find the CDR3b column — it follows CDR-H1/CDR-H2.
  cdr3_col <- intersect(colnames(raw), c("CDR3b", "CDR3", "CDR-H3"))[1]
  if (is.na(cdr3_col)) {
    # Fall back: position 13 per visual header inspection.
    cdr3_col <- colnames(raw)[13]
    message("Using positional CDR3b column: ", cdr3_col)
  }

  raw <- raw %>%
    rename(CDR3b = !!cdr3_col) %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(
      CDR3b,
      TRBV    = .data[["TCRBV"]],
      TRBJ    = .data[["TCRBJ"]],
      patient = .data[["subject"]],
      HLA     = .data[["HLA"]],
      counts  = suppressWarnings(as.integer(.data[["Reads"]]))
    ) %>%
    mutate(counts = ifelse(is.na(counts) | counts < 1, 1L, counts)) %>%
    distinct()

  labels_tbl <- raw %>%
    filter(!is.na(.data[["antigen"]])) %>%
    transmute(
      CDR3b,
      antigen    = .data[["antigen"]],
      peptide    = .data[["peptide"]],
      mhc_allele = .data[["HLA"]]
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

  # MOESM3 sheet "bulk TCR" — clean tabular TCR table.
  # Headers: CDR3b, Vb, Jb, CDR3a, Va, Ja, Individual, Counts
  raw <- read_excel(src, sheet = "bulk TCR", .name_repair = "minimal")

  raw <- raw %>%
    rename(CDR3b = .data[["CDR3b"]]) %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(
      CDR3b,
      TRBV    = .data[["Vb"]],
      TRBJ    = .data[["Jb"]],
      patient = as.character(.data[["Individual"]]),
      HLA     = NA_character_,
      counts  = suppressWarnings(as.integer(.data[["Counts"]]))
    ) %>%
    mutate(counts = ifelse(is.na(counts) | counts < 1, 1L, counts)) %>%
    distinct()

  # Antigen labels for Huang: bulk TCR is Mtb-specific by selection (TB cohort
  # tetramer-sorted), so label all retained rows with antigen = "Mtb". For a
  # finer-grained label, MOESM4 "known_CDR3" lists CDR3+epitope pairs — pull
  # those in for any matching CDR3 to enrich the labels.
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
