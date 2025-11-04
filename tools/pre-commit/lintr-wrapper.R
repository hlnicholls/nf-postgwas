#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
files <- args

is_r_file <- function(f) grepl('\\.(R|r|Rmd|rmd)$', f)
files <- Filter(is_r_file, files)
if (length(files) == 0) {
  quit(status = 0)
}

if (!requireNamespace("lintr", quietly = TRUE)) {
  message("Installing 'lintr' from CRAN (required for pre-commit).")
  install.packages("lintr", repos = "https://cloud.r-project.org")
}

library(lintr)
failed <- FALSE
for (f in files) {
  l <- lintr::lint(f)
  if (length(l) > 0) {
    failed <- TRUE
    cat(paste0("\nLint issues in ", f, "\n"))
    print(l)
  }
}

if (failed) {
  cat('\nR linting failed. Please run lintr locally and fix issues.\n')
  quit(status = 1)
} else {
  quit(status = 0)
}
