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

# Pin immGLIPH at the exact commit under Bioconductor review (v0.99.3 on the
# bioc-review branch). Pinning by SHA keeps results reproducible even if the
# branch moves.
renv::install("github::BorchLab/immGLIPH@c36272c71883259ab3b3d659ab0765194c831ad6")

renv::snapshot(prompt = FALSE)

cat("\nBootstrap complete. renv.lock written.\n")
