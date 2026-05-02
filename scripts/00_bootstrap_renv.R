#!/usr/bin/env Rscript
# One-time R environment bootstrap for the immGLIPH concordance benchmark.
#
# Run from the repo root:
#   Rscript scripts/00_bootstrap_renv.R
#
# After this completes, `renv.lock` will pin every package version. Future
# collaborators only need:  Rscript -e 'renv::restore()'

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::init(bare = TRUE, restart = FALSE)

renv::install(c(
  "tidyverse",
  "mclust",
  "aricode",
  "fs",
  "glue",
  "yaml",
  "jsonlite",
  "optparse",
  "readxl",
  "writexl"
))

renv::install("bioc::BiocParallel")

# Pin immGLIPH at the version under Bioconductor review.
renv::install("github::ncborcherding/immGLIPH@v0.99.3")

renv::snapshot(prompt = FALSE)

cat("\nBootstrap complete. renv.lock written.\n")
