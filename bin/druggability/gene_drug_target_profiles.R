#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(forcats)
  library(R6)
})

# ------------------------------------------------------------------
# Load run-time configuration created by the Nextflow shim
# (defines: prioritised_genes, gene_targets, databases, druggability_output_path, etc.)
# ------------------------------------------------------------------
source("config_R.R")

# ----------------------------- Utilities -----------------------------

safe_fread <- function(path, ...) {
  if (!file.exists(path)) {
    warning(sprintf("⚠️ File not found (skipping): %s", path))
    return(NULL)
  }
  tryCatch(
    data.table::fread(path, ...),
    error = function(e) {
      warning(sprintf("⚠️ Failed to read: %s\nReason: %s", path, e$message))
      NULL
    }
  )
}

safe_merge <- function(x, y, by, all.x = TRUE, allow.cartesian = TRUE) {
  if (is.null(y) || nrow(y) == 0) return(x)
  merge(x, y, by = by, all.x = all.x, allow.cartesian = allow.cartesian)
}

distinct_rows <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  unique(dt)
}

# ----------------------------- R6 Class ------------------------------

GeneDrugTargetProfiles <- R6::R6Class(
  "GeneDrugTargetProfiles",
  public = list(
    # ---- configuration (from config_R.R) ----
    prioritised_genes_path = NULL,
    gene_targets_path      = NULL,
    databases              = NULL,
    druggability_output_path = NULL,

    # optional disease terms (for consistency with other steps)
    disease_terms = NULL,
    disease_regex = NULL,

    initialize = function() {
      # Required paths
      self$prioritised_genes_path   <- get("prioritised_genes",        envir = .GlobalEnv)
      self$gene_targets_path        <- get("gene_targets",             envir = .GlobalEnv)
      self$databases                <- get("databases",                envir = .GlobalEnv)
      self$druggability_output_path <- get("druggability_output_path", envir = .GlobalEnv)

      # Ensure output dir
      dir.create(self$druggability_output_path, showWarnings = FALSE, recursive = TRUE)

      # disease terms/regex from shim (not used in this script’s logic, but we keep consistent warnings)
      default_disease <- "cardio|cardiac|atrial|myocard|arrhyt|vascular|heart|hypertroph|dilated|left ventric|right ventric|hypertension|blood pressure|diabetes|obesity|hypotension|heart failure"

      self$disease_terms <- if (exists("disease_terms", envir = .GlobalEnv)) {
        get("disease_terms", envir = .GlobalEnv)
      } else {
        warning("⚠️ No disease_terms found in config_R.R or params.yaml; using default cardiovascular-related terms.")
        character()
      }

      self$disease_regex <- if (exists("disease_regex", envir = .GlobalEnv)) {
        get("disease_regex", envir = .GlobalEnv)
      } else {
        warning("⚠️ No disease_regex found; using default pattern for cardiovascular-related terms.")
        default_disease
      }
    },

    run = function() {
      message("▶ Loading prioritised genes …")
      prio <- safe_fread(self$prioritised_genes_path)
      if (is.null(prio) || nrow(prio) == 0) {
        stop(sprintf("No data in prioritised genes file: %s", self$prioritised_genes_path))
      }

      message("▶ Loading target classification table …")
      gt <- safe_fread(self$gene_targets_path)
      if (is.null(gt) || nrow(gt) == 0) {
        stop(sprintf("No data in gene target file: %s", self$gene_targets_path))
      }

      # Keep necessary columns from prioritised genes
      prio_keep <- prio %>%
        dplyr::select(
          Nearest_Gene_10kb,
          Gene_Prioritisation_Score,
          `DGIdb_drug(s) (interaction score)`,
          `OT_drug(s)`,
          Druggable_dgidb
        )

      # Filter gene targets to druggable/prioritised gene symbols
      message("▶ Aligning HGNC symbols with prioritised genes …")
      prio_keep <- prio_keep %>% filter(Nearest_Gene_10kb != "" & !is.na(Nearest_Gene_10kb))
      gt_druggable <- gt %>% filter(HGNC_symbol %in% prio_keep$Nearest_Gene_10kb)

      # Match columns and merge
      prio_keep <- prio_keep %>% rename(HGNC_symbol = Nearest_Gene_10kb)
      merged <- safe_merge(gt_druggable, prio_keep, by = "HGNC_symbol")

      # ------------------ External annotations ------------------

      # HIPred
      hipred_path <- file.path(self$databases, "HIPred/HIPred.tsv")
      message("▶ Attaching HIPred (if available) …")
      hipred <- safe_fread(hipred_path)
      merged <- safe_merge(merged, hipred, by = "HGNC_symbol")

      # GDI (keep HGNC_symbol & GDI_Prediction)
      gdi_path <- file.path(self$databases, "gene_damage_index/GDI.txt")
      message("▶ Attaching Gene Damage Index (if available) …")
      gdi <- safe_fread(gdi_path)
      if (!is.null(gdi) && nrow(gdi) > 0) {
        gdi <- gdi %>% dplyr::select(any_of(c("HGNC_symbol", "GDI_Prediction")))
        merged <- safe_merge(merged, gdi, by = "HGNC_symbol")
      }

      # ------------------ Finalize & Write ------------------

      merged <- distinct_rows(merged)

      out_file <- file.path(self$druggability_output_path, "Druggable_gene_targets.csv")
      message(sprintf("💾 Writing: %s", out_file))
      data.table::fwrite(merged, out_file, row.names = FALSE)

      message("✅ Gene-drug target profile generation complete.")
      invisible(TRUE)
    }
  )
)

# ----------------------------- Execute ------------------------------
runner <- GeneDrugTargetProfiles$new()
invisible(runner$run())
