#!/usr/bin/env Rscript

# eQTL colocalisation analysis pipeline (OOP / R6)
# Requirements: coloc, data.table
# Note: GTEx v8 allpairs files are huge; this code uses zcat|awk|sed to read per-gene chunks.
# IMPORTANT: Ensure variant positions are in hg38.
# Inputs expected from config_R.R:
#   - traits: character vector of trait names
#   - gwas_path: path to per-trait GWAS TSVs (<trait>_38_37_rsids.txt)
#   - eqtl_in: directory containing "<trait>_eQTL_candidates.txt"
#   - eqtl_dir: base dir that contains "gtex_eqtl_assocations/<tissue><suffix>"
#   - coloc_path: output directory
#   - all_loci: (kept for compatibility; not required by this step)
# Outputs:
#   - For each trait: "<trait>_eQTL_COLOC.tsv" in coloc_path with PP3/PP4

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(R6)
  library(magrittr)
  library(dplyr)
})

# Source the runroot config file explicitly from the current working directory.
# Nextflow cds into the runroot before executing this script, so load the
# configuration from there. Produce a helpful error if it is missing.
cfg <- file.path(getwd(), "config_R.R")
if (!file.exists(cfg)) {
  stop(sprintf("Missing runroot config: %s\nEnsure CREATE_CONFIG_SHIMS ran and produced config_R.R in the runroot directory.", cfg))
}
source(cfg)

# -----------------------
# User-configurable knobs
# -----------------------
eqtl_cand_dir <- eqtl_in
eqtl_suff     <- ".v8.EUR.allpairs.AllChr.txt.gz"
save_dir      <- coloc_path
lead_var_file <- all_loci
plot_dir      <- coloc_path
bp_lim        <- 500e3         # retained for compatibility (not used in current filtering)
nlines        <- 100e3         # max lines to pull per gene (windowed read)
p12           <- 1e-6          # prior Pr(SNP associated with both traits)
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
setwd(save_dir)

ColocEQTL <- R6::R6Class(
  "ColocEQTL",
  public = list(
    traits = NULL,
    gwas_path = NULL,
    eqtl_cand_dir = NULL,
    eqtl_dir = NULL,
    eqtl_suff = NULL,
    out_dir = NULL,
    nlines = NULL,
    p12 = NULL,

    initialize = function(traits, gwas_path, eqtl_cand_dir, eqtl_dir, eqtl_suff,
                          out_dir, nlines = 100e3, p12 = 1e-6) {
      self$traits         <- traits
      self$gwas_path      <- gwas_path
      self$eqtl_cand_dir  <- eqtl_cand_dir
      self$eqtl_dir       <- eqtl_dir
      self$eqtl_suff      <- eqtl_suff
      self$out_dir        <- out_dir
      self$nlines         <- as.integer(nlines)
      self$p12            <- p12
      dir.create(self$out_dir, showWarnings = FALSE, recursive = TRUE)
    },

    run = function() {
      for (phenotype in self$traits) {
        message("\n========================================")
        message("Processing trait: ", phenotype)
        message("========================================")

        tryCatch(
          {
            gwas <- self$load_gwas(phenotype)
            candidates <- self$load_candidates(phenotype)
            tissues <- unique(candidates$tissue)
            message("Tissues to test: ", paste(tissues, collapse = ", "))

            # Preallocate output columns
            candidates$PP3 <- NA_real_
            candidates$PP4 <- NA_real_

            for (tissue in tissues) {
              message(">> Tissue: ", tissue)
              eqtl_file <- file.path(self$eqtl_dir, "gtex_eqtl_assocations", paste0(tissue, self$eqtl_suff))
              if (!file.exists(eqtl_file)) {
                warning("Missing eQTL file for tissue: ", tissue, " -> ", eqtl_file)
                next
              }
              # Fetch header (column names)
              cols_eqtl <- scan(eqtl_file, what = character(), nlines = 1, sep = "\t", quiet = TRUE)

              idx <- which(candidates$tissue %in% tissue)
              n_done <- 0L

              for (j in idx) {
                gene_id <- candidates$gene_id[j]
                message(sprintf("   [%d/%d] gene_id: %s", n_done + 1L, length(idx), gene_id))

                # Locate first line for this gene in the giant gz (1-based line index in the uncompressed stream)
                line_start <- self$first_line_for_gene(eqtl_file, gene_id)
                if (is.na(line_start)) {
                  warning("   Gene not found in ", basename(eqtl_file), ": ", gene_id)
                  n_done <- n_done + 1L
                  next
                }

                # Read a window (nlines) after the first occurrence
                eqtl_chunk <- self$read_eqtl_chunk(eqtl_file, line_start, self$nlines, cols_eqtl)
                if (is.null(eqtl_chunk) || nrow(eqtl_chunk) == 0) {
                  warning("   Empty eQTL chunk for gene: ", gene_id)
                  n_done <- n_done + 1L
                  next
                }

                # Filter for exact gene
                eqtl_chunk <- subset(eqtl_chunk, complete.cases(eqtl_chunk) & phenotype_id == gene_id)
                if (!nrow(eqtl_chunk)) {
                  warning("   No rows for exact gene_id after filtering: ", gene_id)
                  n_done <- n_done + 1L
                  next
                }

                # Prepare eQTL fields
                eqtl_prep <- self$prepare_eqtl(eqtl_chunk)

                # Align region in GWAS by chr & POS span observed in this eQTL gene block
                tmp_gwas <- self$subset_gwas_for_eqtl_region(gwas, eqtl_prep)
                if (!nrow(tmp_gwas)) {
                  warning("   No GWAS variants overlapping region for gene: ", gene_id)
                  n_done <- n_done + 1L
                  next
                }

                # Merge by exact variant_id (forward and flipped alleles)
                eqtl_for_coloc <- self$merge_eqtl_gwas(tmp_gwas, eqtl_prep)
                if (!nrow(eqtl_for_coloc)) {
                  warning("   No overlapping variants after merging (including flip) for gene: ", gene_id)
                  n_done <- n_done + 1L
                  next
                }

                # Compute coloc.abf
                res <- self$run_coloc(eqtl_for_coloc, tmp_gwas)

                # Store PP3/PP4 (as percentages)
                if (!is.null(res) && !is.null(res$summary)) {
                  candidates$PP3[j] <- round(100 * as.numeric(res$summary[5]))
                  candidates$PP4[j] <- round(100 * as.numeric(res$summary[6]))
                } else {
                  candidates$PP3[j] <- NA_real_
                  candidates$PP4[j] <- NA_real_
                }

                n_done <- n_done + 1L
              } # end per-gene loop
            } # end per-tissue loop

            # Merge locus-positioning info for output and write
            out <- merge(
              candidates,
              gwas[, c("CHROM", "GENPOS", "SNP")],
              by = c("CHROM", "GENPOS"),
              all.x = TRUE
            )

            out_file <- file.path(self$out_dir, paste0(phenotype, "_eQTL_COLOC.tsv"))
            write.table(out, out_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
            message("Wrote: ", out_file)
            message("Successfully completed GTEx coloc for: ", phenotype)
          },
          error = function(e) {
            message("ERROR processing trait ", phenotype, ": ", e$message)
            message("Continuing with next trait...")
          }
        )
      } # end per-trait
    },

    # -----------------------
    # Helpers
    # -----------------------

    load_gwas = function(phenotype) {
      gwas_file <- file.path(self$gwas_path, paste0(phenotype, "_38_37_rsids.txt"))
      if (!file.exists(gwas_file)) stop("GWAS file not found: ", gwas_file)
      gwas <- data.table::fread(gwas_file)

      # Normalise identifiers
      gwas$cptid <- gsub(":", "_", gwas$SNP)

      # Ensure MAF present: derive from A1FREQ if missing
      if (!("MAF" %in% names(gwas))) {
        if ("A1FREQ" %in% names(gwas)) {
          gwas$MAF <- pmin(gwas$A1FREQ, 1 - gwas$A1FREQ)
        } else {
          warning("GWAS does not have MAF or A1FREQ; coloc may fail if MAF is required.")
          gwas$MAF <- NA_real_
        }
      }

      # Ensure numeric N per row if available; otherwise keep as-is
      if ("N" %in% names(gwas)) {
        suppressWarnings(gwas$N <- as.numeric(gwas$N))
      }

      gwas
    },

    load_candidates = function(phenotype) {
      cand_file <- file.path(self$eqtl_cand_dir, paste0("/", phenotype, "_eQTL_candidates.txt"))
      if (!file.exists(cand_file)) stop("Candidate eQTL file not found: ", cand_file)
      read.table(cand_file, header = TRUE, stringsAsFactors = FALSE)
    },

    first_line_for_gene = function(eqtl_file, gene_id) {
      # returns first 1-based line index in the uncompressed stream; NA if not found
      cmd <- sprintf("zcat %s | awk '{print $1}' | grep -w -n -m1 '%s'", shQuote(eqtl_file), gene_id)
      out <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
      if (!length(out)) return(NA_real_)
      as.numeric(strsplit(out, ":", fixed = TRUE)[[1]][1])
    },

    read_eqtl_chunk = function(eqtl_file, start_line, nlines, header_cols) {
      # Read [start_line, start_line + nlines] from gz via sed; skip header line already counted as 1
      cmd <- sprintf("zcat %s | sed -n '%d,%dp'", shQuote(eqtl_file), start_line, start_line + nlines)
      df <- tryCatch(
        data.table::fread(cmd, header = FALSE, sep = "\t", showProgress = FALSE),
        warning = function(w) suppressWarnings(data.table::fread(cmd, header = FALSE, sep = "\t", showProgress = FALSE)),
        error = function(e) NULL
      )
      if (is.null(df) || !nrow(df)) return(NULL)
      colnames(df) <- header_cols
      df
    },

    prepare_eqtl = function(eqtl_chunk) {
      # Compute N_eQTL and allele fields; ensure numeric types
      eq <- as.data.frame(eqtl_chunk)
      # N_eQTL: make a pragmatic estimate; original used round(ma_samples * 1/maf)
      if (!("ma_samples" %in% names(eq)) || !("maf" %in% names(eq))) {
        stop("Required GTEx columns missing in eQTL: ma_samples or maf")
      }
      eq$N_eQTL <- round(eq$ma_samples * (1 / eq$maf))
      eq$ma_samples <- NULL
      colnames(eq)[colnames(eq) == "maf"] <- "MAF_eQTL"

      # Parse variant_id "chrCHR_POS_REF_ALT_b38"
      parts <- strsplit(eq$variant_id, "_", fixed = TRUE)
      eq$CHR <- as.numeric(gsub("chr", "", vapply(parts, `[`, character(1), 1L)))
      eq$POS <- as.numeric(vapply(parts, `[`, character(1), 2L))
      eq$REF <- vapply(parts, `[`, character(1), 3L)
      eq$ALT <- vapply(parts, `[`, character(1), 4L)

      eq
    },

    subset_gwas_for_eqtl_region = function(gwas, eq) {
      chr <- unique(eq$CHR)
      if (length(chr) != 1L || is.na(chr)) return(gwas[0, ])
      pos_min <- min(eq$POS, na.rm = TRUE)
      pos_max <- max(eq$POS, na.rm = TRUE)
      subset(gwas, CHROM == chr & GENPOS >= pos_min & GENPOS <= pos_max)
    },

    merge_eqtl_gwas = function(tmp_gwas, eq) {
      # Forward orientation variant_id (as in GTEx)
      gwas_fwd <- tmp_gwas
      gwas_fwd$variant_id <- paste0(
        "chr",
        paste(gwas_fwd$CHROM, gwas_fwd$GENPOS, gwas_fwd$ALLELE1, gwas_fwd$ALLELE0, "b38", sep = "_")
      )

      # Merge forward
      m_fwd <- merge(
        gwas_fwd,
        eq[, c("variant_id", "MAF_eQTL", "slope", "slope_se", "N_eQTL", "pval_nominal")],
        by = "variant_id"
      )

      # Flipped orientation (swap REF/ALT)
      eq_flip <- eq
      eq_flip$variant_id <- paste0("chr", paste(eq$CHR, eq$POS, eq$ALT, eq$REF, "b38", sep = "_"))

      m_flip <- merge(
        gwas_fwd,
        eq_flip[, c("variant_id", "MAF_eQTL", "slope", "slope_se", "N_eQTL", "pval_nominal")],
        by = "variant_id"
      )

      # Combine and order
      eqtl_coloc <- rbind(m_fwd, m_flip)
      if (!nrow(eqtl_coloc)) return(eqtl_coloc)
      eqtl_coloc <- eqtl_coloc[order(eqtl_coloc$CHROM, eqtl_coloc$GENPOS), ]
      eqtl_coloc
    },

    run_coloc = function(eqtl_coloc, tmp_gwas) {
      # Choose N for GWAS—use max observed N in the overlapping window, else first row fallback
      N_gwas <- if ("N" %in% names(tmp_gwas) && any(!is.na(tmp_gwas$N))) {
        max(tmp_gwas$N, na.rm = TRUE)
      } else {
        tmp_gwas$N[1]
      }

      # Build coloc datasets. sdY=1 for standardised quantitative trait summary stats.
      # GWAS MAF: try from tmp_gwas if present after merge; if not, it's acceptable to pass NULL.
      maf_gwas <- if ("MAF" %in% names(eqtl_coloc)) eqtl_coloc$MAF else NULL

      tryCatch(
        coloc.abf(
          dataset1 = list(
            sdY    = 1,
            beta   = eqtl_coloc$BETA,
            varbeta= eqtl_coloc$SE^2,
            N      = N_gwas,
            type   = "quant",
            MAF    = maf_gwas,
            snp    = eqtl_coloc$SNP
          ),
          dataset2 = list(
            beta    = eqtl_coloc$slope,
            varbeta = eqtl_coloc$slope_se^2,
            N       = max(eqtl_coloc$N_eQTL, na.rm = TRUE),
            type    = "quant",
            MAF     = eqtl_coloc$MAF_eQTL,
            snp     = eqtl_coloc$SNP
          ),
          p12 = self$p12
        ),
        error = function(e) {
          warning("coloc.abf failed: ", e$message)
          NULL
        }
      )
    }
  )
)

# ----------------------------- Execute ------------------------------
pipeline <- ColocEQTL$new(
  traits         = traits,
  gwas_path      = gwas_path,
  eqtl_cand_dir  = eqtl_cand_dir,
  eqtl_dir       = eqtl_dir,
  eqtl_suff      = eqtl_suff,
  out_dir        = save_dir,
  nlines         = nlines,
  p12            = p12
)
pipeline$run()
