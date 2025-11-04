#!/usr/bin/env Rscript

# Prepare GWAS summary statistics for MAGMA analysis (single trait)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(here)
})

# Source shared config
source("config_R.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: prep_SNPs_for_magma.R <trait> <is_mtag>")
}

trait <- args[1]
is_mtag <- args[2]

# Optional: third arg can be the path to the input GWAS file (rsids file)
if (length(args) >= 3) {
  full_path <- args[3]
} else {
  # Fallback: construct expected filename from trait and gwas_input_path
  file_name <- paste0(trait, "_38_37_rsids.txt")
  full_path <- file.path(gwas_input_path, file_name)
}

# Output filename expected by downstream MAGMA wrapper
file_name_out <- paste0(trait, "_38_rsid_magma.txt")

if (!file.exists(full_path)) {
  stop(sprintf("Input file not found: %s", full_path))
}

gwas <- fread(full_path)

gwas <- gwas %>%
  select(rsid_1kg, SNP, P, CHROM, GENPOS, ALLELE0, ALLELE1, N, MAF, BETA, SE) %>%
  filter(!is.na(rsid_1kg) & rsid_1kg != "")

colnames(gwas)[c(1, 2)] <- c('SNP', 'SNP2')

full_path_out <- file.path(magma_output_path, file_name_out)
cat(sprintf("Writing: %s (%d SNPs)\n", full_path_out, nrow(gwas)))
fwrite(gwas, full_path_out, sep = "\t")
