#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(glue)
  library(here)
})

source("config_R.R")

# ------------------------------------------------------------------------------------
# Identify lead snps with close p-values in other GWAS summary stats
# ------------------------------------------------------------------------------------
GWASClosePChecker <- R6::R6Class(
  "GWASClosePChecker",
  public = list(
    all_loci      = NULL,
    traits        = NULL,
    gwas_path     = NULL,
    threshold     = 5e-5,
    output_table  = NULL,

    initialize = function(all_loci, traits, gwas_path, threshold = 5e-5) {
      # Store config values
      self$all_loci  <- all_loci
      self$traits    <- traits
      self$gwas_path <- gwas_path
      self$threshold <- threshold
    },

    load_loci = function() {
      # Load loci file and coerce P to numeric, create Has_Close_P as NA
      self$output_table <- fread(self$all_loci)
      self$output_table$P <- as.numeric(self$output_table$P)
      self$output_table$Has_Close_P <- NA
      invisible(self)
    },

    process_trait = function(phenotype) {
      # Load GWAS summary stats for the phenotype and filter by threshold
      gwas_file <- glue("{self$gwas_path}/{phenotype}_38_37.txt")
      gwas <- fread(gwas_file)
      gwas <- dplyr::filter(gwas, P < self$threshold)

      # Preserve original behavior: print the phenotype name
      print(phenotype)

      # For each row in output_table (same loop style and checks as original)
      for (i in seq_len(nrow(self$output_table))) {
        if (self$output_table$Phenotype[i] != phenotype) {
          matching_snp <- dplyr::filter(gwas, SNP == self$output_table$SNP[i])
          if (nrow(matching_snp) > 0) {
            p_value <- formatC(min(matching_snp$P), format = "e", digits = 3)
            if (is.na(self$output_table$Has_Close_P[i])) {
              self$output_table$Has_Close_P[i] <- glue("Yes({p_value},{phenotype}, Single-trait)")
            } else {
            }
          }
        }
      }
      invisible(self)
    },

    run = function() {
      # Load loci table
      self$load_loci()

      # Loop through all phenotypes (traits)
      for (phenotype in self$traits) {
        self$process_trait(phenotype)
      }

      # Write the updated loci file back to the same location
      fwrite(self$output_table, self$all_loci)
      invisible(self)
    }
  )
)

# ----------------------------- Execute -----------------------------
GWASClosePChecker$
  new(all_loci = all_loci,
      traits   = traits,
      gwas_path = gwas_path,
      threshold = 5e-5)$
  run()
