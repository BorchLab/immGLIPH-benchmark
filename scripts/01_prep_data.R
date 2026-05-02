#!/usr/bin/env Rscript
# Normalize Glanville 2017 and Huang 2020 supplementary data to the common
# input schema consumed by every downstream tool (immGLIPH, GLIPH Perl, GLIPH2
# binary).
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
  # The Glanville 2017 supplement has multiple sheets across two XLSX files.
  # The TCR table (Supp Table 1) has columns roughly:
  #   CDR3b, TRBV, TRBJ, Subject, Antigen, HLA, Count
  # When the raw file is in place, replace the path glob below.
  xlsx_files <- dir_ls(raw_dir, glob = "*.xlsx")
  if (length(xlsx_files) == 0) {
    stop("No XLSX in ", raw_dir,
         ". Download Glanville 2017 supplementary tables first (see SOURCE.md).")
  }

  # TODO: once raw file is in place, set the correct file + sheet here.
  tcr_path <- xlsx_files[1]
  raw <- read_excel(tcr_path, sheet = 1)

  # Column normalization. Adjust on first run if Glanville's actual headers differ.
  col_map <- c(
    CDR3b   = "CDR3b",
    TRBV    = "TRBV",
    TRBJ    = "TRBJ",
    patient = "Subject",
    HLA     = "HLA",
    counts  = "Count",
    antigen = "Antigen",
    peptide = "Epitope",
    mhc_allele = "MHC"
  )
  missing <- setdiff(col_map, colnames(raw))
  if (length(missing) > 0) {
    warning("Columns not found, will be NA: ",
            paste(missing, collapse = ", "))
  }
  raw <- raw %>%
    rename(any_of(col_map)) %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(CDR3b,
              TRBV   = if ("TRBV" %in% names(.)) TRBV else NA_character_,
              TRBJ   = if ("TRBJ" %in% names(.)) TRBJ else NA_character_,
              patient = if ("patient" %in% names(.)) patient else NA_character_,
              HLA    = if ("HLA" %in% names(.)) HLA else NA_character_,
              counts = if ("counts" %in% names(.)) counts else 1L) %>%
    distinct()

  labels_tbl <- raw %>%
    filter(!is.na(antigen)) %>%
    transmute(CDR3b, antigen,
              peptide = if ("peptide" %in% names(.)) peptide else NA_character_,
              mhc_allele = if ("mhc_allele" %in% names(.)) mhc_allele else NA_character_) %>%
    distinct()

  out_dir <- "data/glanville2017"
  write_tsv(input_tbl, file.path(out_dir, "cdr3_input.tsv"))
  write_tsv(labels_tbl, file.path(out_dir, "antigen_labels.tsv"))

  message(sprintf("Glanville: %d input sequences, %d labeled.",
                  nrow(input_tbl), nrow(labels_tbl)))
}

prep_huang <- function() {
  raw_dir <- "data/huang2020/raw"
  xlsx_files <- dir_ls(raw_dir, glob = "*.xlsx")
  if (length(xlsx_files) == 0) {
    stop("No XLSX in ", raw_dir,
         ". Download Huang 2020 supplementary tables first (see SOURCE.md).")
  }

  # TODO: once raw file is in place, set the correct file + sheet here.
  # Huang 2020 TCR table column names (typically): cdr3, V, J, subject, HLA,
  # frequency, antigen.species, antigen.epitope, mhc.a
  tcr_path <- xlsx_files[1]
  raw <- read_excel(tcr_path, sheet = 1)

  col_map <- c(
    CDR3b   = "cdr3",
    TRBV    = "V",
    TRBJ    = "J",
    patient = "subject",
    HLA     = "HLA",
    counts  = "frequency",
    antigen = "antigen.species",
    peptide = "antigen.epitope",
    mhc_allele = "mhc.a"
  )
  missing <- setdiff(col_map, colnames(raw))
  if (length(missing) > 0) {
    warning("Columns not found, will be NA: ",
            paste(missing, collapse = ", "))
  }
  raw <- raw %>%
    rename(any_of(col_map)) %>%
    filter(valid_cdr3(CDR3b))

  input_tbl <- raw %>%
    transmute(CDR3b,
              TRBV   = if ("TRBV" %in% names(.)) TRBV else NA_character_,
              TRBJ   = if ("TRBJ" %in% names(.)) TRBJ else NA_character_,
              patient = if ("patient" %in% names(.)) patient else NA_character_,
              HLA    = if ("HLA" %in% names(.)) HLA else NA_character_,
              counts = if ("counts" %in% names(.)) counts else 1L) %>%
    distinct()

  labels_tbl <- raw %>%
    filter(!is.na(antigen)) %>%
    transmute(CDR3b, antigen,
              peptide = if ("peptide" %in% names(.)) peptide else NA_character_,
              mhc_allele = if ("mhc_allele" %in% names(.)) mhc_allele else NA_character_) %>%
    distinct()

  out_dir <- "data/huang2020"
  write_tsv(input_tbl, file.path(out_dir, "cdr3_input.tsv"))
  write_tsv(labels_tbl, file.path(out_dir, "antigen_labels.tsv"))

  message(sprintf("Huang: %d input sequences, %d labeled.",
                  nrow(input_tbl), nrow(labels_tbl)))
}

switch(opt$dataset,
  glanville = prep_glanville(),
  huang     = prep_huang()
)
