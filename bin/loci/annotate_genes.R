#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(vautils)
  library(here)
})

source("config_R.R")

# Annotates loci (lead + proxies) with nearest genes within 10kb
LDGeneAnnotator <- R6::R6Class(
  classname = "LDGeneAnnotator",
  public = list(
    # ---- fields ----
    output_path   = NULL,
    gwas_file     = NULL,
    output_file   = NULL,
    build         = NULL,
    flank_kb      = NULL,
    col_snp       = NULL,
    col_chr       = NULL,
    col_bp        = NULL,

    # ---- ctor ----
    initialize = function(output_path,
                          build      = "hg38",
                          flank_kb   = 10,
                          col_snp    = "SNP",
                          col_chr    = "CHROM",
                          col_bp     = "GENPOS") {

      stopifnot(!is.null(output_path), nzchar(output_path))

      self$output_path <- normalizePath(output_path, mustWork = FALSE)
  self$gwas_file   <- file.path(self$output_path, "Loci_Preprocessing", "all_traits_loci_38_with_ld.txt")
  self$output_file <- file.path(self$output_path, "Loci_Preprocessing", "all_traits_loci_38_with_ld_genes.txt")

      self$build   <- build
      self$flank_kb <- flank_kb
      self$col_snp <- col_snp
      self$col_chr <- col_chr
      self$col_bp  <- col_bp
    },

    # ---- I/O helpers ----
    read_gwas = function() {
      if (!file.exists(self$gwas_file)) {
        stop(sprintf("Input GWAS file not found: %s", self$gwas_file))
      }
      data.table::fread(self$gwas_file)
    },

    ensure_outdir = function() {
      outdir <- dirname(self$output_file)
      if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    },

    write_output = function(dt) {
      self$ensure_outdir()
      data.table::fwrite(dt, file = self$output_file, sep = "\t", row.names = FALSE)
    },

    # ---- core annotate ----
    annotate_genes = function(gwas_dt) {
      # find nearest gene within ±flank_kb kilobases
      geneN <- vautils::find_nearest_gene(
        gwas_dt,
        flanking = self$flank_kb,
        build    = self$build,
        collapse = FALSE,
        snp      = self$col_snp,
        chr      = self$col_chr,
        bp       = self$col_bp
      )

      geneN <- geneN %>% dplyr::rename(SNP = .data$rsid)
      geneN <- geneN %>% dplyr::mutate(distance = ifelse(.data$distance == "intergenic", 0, .data$distance))

      # one gene per SNP: pick row with smallest |distance|
      gene_unique <-
        geneN[order(geneN$SNP, abs(as.numeric(geneN$distance))), ] %>%
        subset(!duplicated(SNP)) %>%
        dplyr::select(SNP, GENE)

      # merge and put GENE first, named Gene_Symbol, then everything else
      merged <-
        merge(gene_unique, gwas_dt, by = "SNP") %>%
        dplyr::select(Gene_Symbol = .data$GENE, dplyr::everything())

      merged
    },

    # ---- orchestrator ----
    run = function() {
      message("Getting genes within 10kb for all loci (leads and proxies) for all single traits...")
      gwas_dt <- self$read_gwas()
      annotated <- self$annotate_genes(gwas_dt)
      self$write_output(annotated)
      message(sprintf("Wrote: %s", self$output_file))
      invisible(self$output_file)
    }
  )
)

# ----------------------------- Execute -----------------------------
annotator <- LDGeneAnnotator$new(
  output_path = output_path,
  build       = "hg38",
  flank_kb    = 10,
  col_snp     = "SNP",
  col_chr     = "CHROM",
  col_bp      = "GENPOS"
)

annotator$run()
