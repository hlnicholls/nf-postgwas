#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(R6)
})

# ------------------------------------------------------------
# Load shared config written by the Nextflow shim into runroot
# (must define: databases, var_file, var_path; optionally HiC_tissues)
# ------------------------------------------------------------
source("config_R.R")

HiCAnnotator <- R6::R6Class(
  "HiCAnnotator",
  public = list(
    # config
    databases = NULL,
    var_file  = NULL,
    var_path  = NULL,

    # inputs
    hic_tissues = NULL,
    chrom_file  = NULL,   # PCHi-C file path (hg19)

    # ctor
    initialize = function() {
      # Required paths from config
      self$databases <- get("databases", envir = .GlobalEnv)
      self$var_file  <- get("var_file",  envir = .GlobalEnv)
      self$var_path  <- get("var_path",  envir = .GlobalEnv)

      # Hi-C tissues (optional param set in shim as HiC_tissues)
      default_tissues <- c("Right Atrium", "Right Ventricle", "Left Ventricle", "Aorta")
      if (exists("HiC_tissues", envir = .GlobalEnv)) {
        self$hic_tissues <- get("HiC_tissues", envir = .GlobalEnv)
        if (length(self$hic_tissues) == 0) {
          warning("HiC_tissues is defined but empty; defaulting to cardiac/aortic tissues.")
          self$hic_tissues <- default_tissues
        }
      } else {
        warning("No HiC_tissues found in config_R.R; defaulting to cardiac/aortic tissues.")
        self$hic_tissues <- default_tissues
      }

      # Promoter capture Hi-C (hg19) file path (Jung et al.)
      self$chrom_file <- file.path(self$databases, "hic/41588_2019_494_MOESM3_ESM.txt.gz")
    },

    run = function() {
      # 1) Load PCHi-C and filter by tissues
      pcHiC <- self$load_pcHiC_()
      message("Available tissue counts in PCHi-C (post-filter):")
      print(table(pcHiC$Tissue_type))

      # 2) Build GRanges for interacting fragments (hg19)
      intGR <- self$build_interaction_gr_(pcHiC)

      # 3) Load GWAS/LD table (must contain hg19 coords & alleles)
      ld <- self$load_ld_table_()

      # 4) Overlap SNPs with PCHi-C fragments (hg19)
      anno <- self$annotate_snps_(ld, intGR)

      # 5) Rebuild SNP ID in GRCh38 form and finish
      out <- self$finalize_output_(ld, anno)

      # 6) Write
      out_file <- file.path(self$var_path, "all_traits_in_ld_HiC.tsv")
      fwrite(out, out_file, sep = "\t", quote = FALSE)
      message("Annotated with HiC data → ", out_file)
    },

    # ---------------- internal helpers ----------------

    load_pcHiC_ = function() {
      if (!file.exists(self$chrom_file)) {
        stop("PCHi-C file not found: ", self$chrom_file)
      }
      po <- data.table::fread(self$chrom_file, sep = "\t", header = TRUE, skip = 1)
      # Filter tissues to the configured set
      po %>% dplyr::filter(.data$Tissue_type %in% self$hic_tissues)
    },

    build_interaction_gr_ = function(ints) {
      # Expect columns Interacting_fragment like "chr.start.end" and Promoter, Tissue_type
      req_cols <- c("Interacting_fragment", "Promoter", "Tissue_type")
      missing <- setdiff(req_cols, names(ints))
      if (length(missing) > 0) {
        stop("PCHi-C file is missing required columns: ", paste(missing, collapse = ", "))
      }

      # Parse "chr.start.end"
      parts <- strsplit(ints$Interacting_fragment, "\\.")
      seqnames <- vapply(parts, function(x) x[[1]], "", USE.NAMES = FALSE)
      starts   <- suppressWarnings(as.numeric(vapply(parts, function(x) x[[2]], "", USE.NAMES = FALSE)))
      ends     <- suppressWarnings(as.numeric(vapply(parts, function(x) x[[3]], "", USE.NAMES = FALSE)))

      if (anyNA(starts) || anyNA(ends)) {
        stop("Failed to parse Interacting_fragment coordinates.")
      }

      GRanges(
        seqnames = seqnames,
        ranges   = IRanges(start = starts, end = ends),
        promoter = ints$Promoter,
        tissue   = ints$Tissue_type
      )
    },

    load_ld_table_ = function() {
      if (!file.exists(self$var_file)) {
        stop("var_file not found: ", self$var_file)
      }
      dt <- data.table::fread(self$var_file)

      # Require hg19 coordinates to overlap with hg19 PCHi-C
      req_cols <- c("CHROM", "GENPOS_hg19", "ALLELE0", "ALLELE1")
      missing  <- setdiff(req_cols, names(dt))
      if (length(missing) > 0) {
        stop("Input loci file is missing required columns: ", paste(missing, collapse = ", "),
             ". Ensure liftover to hg19 is present (GENPOS_hg19).")
      }

      # Keep only rows with valid hg19 positions
      dt <- dt %>% filter(!is.na(.data$CHROM) & !is.na(.data$GENPOS_hg19))

      # Normalize chromosome with "chr" prefix for GRanges
      chr_clean <- sub("^chr", "", as.character(dt$CHROM))
      dt$chr_hg19 <- paste0("chr", chr_clean)

      dt
    },

    annotate_snps_ = function(ld, intGR) {
      # Build hg19 SNP GRanges using GENPOS_hg19
      snpGR <- GRanges(
        seqnames = ld$chr_hg19,
        ranges   = IRanges(start = ld$GENPOS_hg19, end = ld$GENPOS_hg19)
      )

      hits <- findOverlaps(snpGR, intGR)
      if (length(hits) == 0) {
        message("No overlaps found between SNPs and PCHi-C fragments with selected tissues.")
        return(data.frame(SNP = character(), HiC_gene = character(), HiC_tissue = character()))
      }

      # Extract annotations
      df <- data.frame(
        idx = as.integer(queryHits(hits)),
        gene = mcols(intGR)$promoter[subjectHits(hits)],
        tissue = mcols(intGR)$tissue[subjectHits(hits)],
        stringsAsFactors = FALSE
      )

      # Build hg19-style SNP ID for merging: CHROM:GENPOS_hg19:ALLELE0:ALLELE1
      # Note: the final output will be re-keyed to GRCh38 SNP IDs later
      df$SNP_hg19 <- paste0(
        sub("^chr", "", ld$chr_hg19[df$idx]), ":",
        ld$GENPOS_hg19[df$idx], ":",
        ld$ALLELE0[df$idx], ":",
        ld$ALLELE1[df$idx]
      )

      out <- df %>%
        dplyr::select(SNP_hg19, HiC_gene = gene, HiC_tissue = tissue) %>%
        distinct()

      out
    },

    finalize_output_ = function(ld, anno_hg19) {
      # If we have no annotations, return empty table with expected columns
      if (nrow(anno_hg19) == 0) {
        return(data.frame(SNP = character(), HiC_gene = character(), HiC_tissue = character()))
      }

      # Build matching hg19 SNP id in LD table to merge with anno_hg19
      ld$SNP_hg19 <- paste0(
        sub("^chr", "", ld$chr_hg19), ":",
        ld$GENPOS_hg19, ":",
        ld$ALLELE0, ":",
        ld$ALLELE1
      )

      # We also need GRCh38 SNP ID in the final file: CHROM:GENPOS:ALLELE0:ALLELE1
      # Require GENPOS (hg38) to exist
      if (!("GENPOS" %in% names(ld))) {
        stop("Input loci file is missing GENPOS (GRCh38 coordinate) required to build final SNP IDs.")
      }

      # Build GRCh38 SNP ID and strip any "chr" prefix
      ld$SNP_hg38 <- paste0(
        sub("^chr", "", as.character(ld$CHROM)), ":",
        ld$GENPOS, ":",
        ld$ALLELE0, ":",
        ld$ALLELE1
      )

      merged <- merge(
        anno_hg19,
        ld[, c("SNP_hg19", "SNP_hg38")],
        by = "SNP_hg19",
        all.x = TRUE
      )

      merged <- merged %>%
        dplyr::select(SNP = SNP_hg38, HiC_gene, HiC_tissue) %>%
        distinct()

      # Clean any residual "chr" prefix that might slip into SNP (defensive)
      merged$SNP <- sub("^chr", "", merged$SNP)

      merged
    }
  )
)

# ----------------------------- Execute -----------------------------
HiCAnnotator$new()$run()
