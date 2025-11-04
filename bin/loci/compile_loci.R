#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(vautils)
  library(data.table)
  library(here)
  library(R6)
})

source("config_R.R")

# ---------- helpers ----------
get_if_exists <- function(name, env = parent.frame()) {
  if (exists(name, envir = env, inherits = TRUE)) get(name, envir = env, inherits = TRUE) else NULL
}

# Prefer gwas_input_path if present; else use gwas_path
.gwas_input_path <- get_if_exists("gwas_input_path", environment())
.gwas_path       <- get_if_exists("gwas_path",       environment())
if (is.null(.gwas_path) && is.null(.gwas_input_path)) {
  stop("Neither gwas_input_path nor gwas_path is defined in config_R.R")
}
resolved_gwas_path <- if (!is.null(.gwas_input_path)) .gwas_input_path else .gwas_path

# Ensure output directory exists
dir.create(file.path(output_path, "Loci_Preprocessing"), showWarnings = FALSE, recursive = TRUE)

# ---------- distance-based pruning function ----------
function.doPruning <- function(dat.unpruned, col_pval, col_CHR, col_BP, BPlim){
  dat.unpruned <- dat.unpruned[order(dat.unpruned[, col_CHR], dat.unpruned[, col_BP]), ]
  dat.pruned   <- dat.unpruned[0, ]
  n <- 0
  for (chr in unique(dat.unpruned[, col_CHR])) {
    dat.unpruned.chr <- subset(dat.unpruned, dat.unpruned[, col_CHR] == chr)
    dat.unpruned.chr$locus <- 0
    dat.unpruned.chr$locus[which(diff(dat.unpruned.chr[, col_BP]) > BPlim) + 1] <- 1
    dat.unpruned.chr$locus <- cumsum(dat.unpruned.chr$locus)
    dat.unpruned.chr <- dat.unpruned.chr[order(dat.unpruned.chr[, col_pval]), ]
    dat.unpruned.chr$locus <- dat.unpruned.chr$locus + n + 1
    dat.pruned <- rbind(dat.pruned, dat.unpruned.chr[!duplicated(dat.unpruned.chr$locus), ])
    n <- max(dat.pruned$locus)
  }
  dat.pruned <- dat.pruned[order(dat.pruned$locus), ]
  dat.pruned
}

# ---------- Compile Loci ----------
LociCompiler <- R6::R6Class(
  "LociCompiler",
  public = list(
    traits = NULL,
    gwas_path = NULL,
    output_path = NULL,

    per_trait_results = NULL,
    combined_loci_df = NULL,

    initialize = function(traits, gwas_path, output_path) {
      self$traits      <- traits
      self$gwas_path   <- gwas_path
      self$output_path <- output_path
    },

    compile_one = function(phenotype, suffix = "_38_37_rsids.txt") {
      message("Processing phenotype: ", phenotype)

      # Read GWAS
      gwas_file <- file.path(self$gwas_path, paste0(phenotype, suffix))
      if (!file.exists(gwas_file))
        stop("Missing GWAS file for trait ", phenotype, ": ", gwas_file)
      test_regenie <- data.table::fread(gwas_file)

      # Filter genome-wide significant
      regenie_filtered_gws <- dplyr::filter(test_regenie, P < 5e-8)
      if (nrow(regenie_filtered_gws) == 0) {
        warning("Phenotype ", phenotype, " has no SNPs with P < 5e-8. Skipping.")
        return(NULL)
      }

      # Prune (1 Mb window)
      pruned <- function.doPruning(
        as.data.frame(regenie_filtered_gws),
        col_pval = "P",
        col_CHR  = "CHROM",
        col_BP   = "GENPOS",
        BPlim    = 1e6
      )

      # Nearest gene (hg38)
      geneN <- find_nearest_gene(
        pruned, flanking = 1000, build = "hg38",
        collapse = FALSE, snp = "SNP", chr = "CHROM", bp = "GENPOS"
      ) %>%
        dplyr::rename(SNP = rsid) %>%
        dplyr::mutate(distance = ifelse(distance == "intergenic", 0, distance))

      # Merge + tidy
      pruned2 <- merge(
        subset(geneN[order(geneN$SNP, abs(as.numeric(geneN$distance))), ], !duplicated(SNP))[, c("SNP", "GENE")],
        pruned,
        by = "SNP"
      ) %>%
        dplyr::arrange(locus) %>%
        dplyr::select(Locus_n = locus, Locus_name = GENE, dplyr::everything())

      cols_to_select <- c(
        "Locus_n", "Locus_name", "SNP", "CHROM", "GENPOS",
        "ALLELE1", "ALLELE0", "A1FREQ", "MAF", "BETA", "SE",
        "P", "GENPOS_hg19", "rsid_1kg"
      )
      pruned2 %<>% dplyr::select(dplyr::all_of(cols_to_select))
      pruned2
    },

    compile_all = function() {
      message("Compiling loci across all phenotypes...")
      results <- purrr::map(self$traits, ~ self$compile_one(.x))
      # keep non-null and name by trait
      keep_idx <- !purrr::map_lgl(results, is.null)
      self$per_trait_results <- results[keep_idx]
      names(self$per_trait_results) <- self$traits[keep_idx]

      if (length(self$per_trait_results) == 0) {
        warning("No traits produced significant loci. Nothing to write.")
        self$combined_loci_df <- tibble()
        return(invisible(self))
      }

      self$combined_loci_df <- dplyr::bind_rows(self$per_trait_results, .id = "Phenotype")

      # Original behavior: rename 4th column to "SNP" (defensive)
      if (ncol(self$combined_loci_df) >= 4) {
        colnames(self$combined_loci_df)[4] <- "SNP"
      }

      self$combined_loci_df$Method <- "Single-trait"
      invisible(self)
    },

    save_outputs = function() {
  outdir <- file.path(self$output_path, "Loci_Preprocessing")
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

      # Full combined
      readr::write_csv(self$combined_loci_df, file.path(outdir, "Singletrait_all_loci.csv"))

      # Unique by SNP
      unique_df <- self$combined_loci_df %>% dplyr::distinct(SNP, .keep_all = TRUE)
      readr::write_csv(unique_df, file.path(outdir, "Singletrait_loci_unique_across_all_traits.csv"))

      if ("SNP" %in% names(self$combined_loci_df)) {
        message("Number of unique Lead SNPs: ", dplyr::n_distinct(self$combined_loci_df$SNP))
      }
      invisible(self)
    },

    run = function() {
      self$compile_all()$save_outputs()
      invisible(self)
    }
  )
)

# ----------------------------- Execute -----------------------------
LociCompiler$
  new(
    traits      = traits,
    gwas_path   = resolved_gwas_path,
    output_path = output_path
  )$
  run()
