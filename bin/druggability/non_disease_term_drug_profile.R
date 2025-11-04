#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(R6)
})

# ------------------------------------------------------------
# Load run-time config (created by the Nextflow shim)
# Defines: databases, druggability_output_path, disease_terms, disease_regex, etc.
# ------------------------------------------------------------
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
  if (is.null(x) || nrow(x) == 0) return(x)
  if (is.null(y) || nrow(y) == 0) return(x)
  merge(x, y, by = by, all.x = all.x, allow.cartesian = allow.cartesian)
}

# ----------------------------- R6 Class ------------------------------

NonDiseaseTermDrugProfile <- R6::R6Class(
  "NonDiseaseTermDrugProfile",
  public = list(
    databases = NULL,
    druggability_output_path = NULL,

    # disease terms / regex
    disease_terms = NULL,
    disease_regex = NULL,

    # input / output files
    drug_table_file = NULL,
    ot_safety_file  = NULL,

    out_hf_file                 = NULL,  # keep 'disease' in filename as requested
    out_ndt_all_file            = NULL,  # nondiseaseterm all
    out_ndt_cvdonly_file        = NULL,  # nondiseaseterm subset where FDA AEs match disease terms

    initialize = function() {
      # Required config
      self$databases                <- get("databases",                envir = .GlobalEnv)
      self$druggability_output_path <- get("druggability_output_path", envir = .GlobalEnv)

      dir.create(self$druggability_output_path, showWarnings = FALSE, recursive = TRUE)

      # Inputs
      self$drug_table_file <- file.path(self$druggability_output_path, "Druggability_table.csv")
      self$ot_safety_file  <- file.path(self$databases, "opentargets/Drug_target_safety_profile.csv")

      # Outputs (names per your request / modules)
      self$out_hf_file          <- file.path(self$druggability_output_path, "disease_drug_safety_FDA_adverse_effects.csv")
      self$out_ndt_all_file     <- file.path(self$druggability_output_path, "nondiseaseterm_gene-drug_OT_FDA_adverse_effects.csv")
      self$out_ndt_cvdonly_file <- file.path(self$druggability_output_path, "nondiseaseterm_gene-drug_OT_FDA_adverse_effects.csv")

      # Disease regex from shim (fallback with warnings)
      default_disease <- "cardio|cardiac|atrial|myocard|arrhyt|vascular|heart|hypertroph|dilated|left ventric|right ventric|hypertension|hypotension|blood pressure|diabetes|obesity|electro|heart failure|cardiac failure"

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
      message("▶ Loading Druggability_table.csv …")
      drug_df <- safe_fread(self$drug_table_file)
      if (is.null(drug_df) || nrow(drug_df) == 0) {
        stop(sprintf("No data in druggability table: %s", self$drug_table_file))
      }

      # maintain compatibility: some runs may or may not have toxicityClass column
      if ("toxicityClass" %in% names(drug_df)) {
        drug_df <- dplyr::select(drug_df, -toxicityClass)
      }

      # column naming consistency from drug table step
      if (!("Interacts_with_non_disease_terms_drug_only" %in% names(drug_df))) {
        # Backwards compatibility with older name (singular 'term')
        old_col <- "Interacts_with_non_disease_term_drug_only"
        if (old_col %in% names(drug_df)) {
          drug_df <- dplyr::rename(drug_df, Interacts_with_non_disease_terms_drug_only = !!old_col)
        } else {
          warning("⚠️ Column 'Interacts_with_non_disease_terms_drug_only' not found; creating 'No' defaults.")
          drug_df$Interacts_with_non_disease_terms_drug_only <- "No"
        }
      }

      message("▶ Loading OpenTargets drug safety profile …")
      ot_safety <- safe_fread(self$ot_safety_file)
      if (is.null(ot_safety) || nrow(ot_safety) == 0) {
        warning(sprintf("⚠️ No safety data available in: %s. Outputs will be empty.", self$ot_safety_file))
        # Write empty shell outputs to keep Nextflow happy
        data.table::fwrite(data.frame(), self$out_hf_file)
        data.table::fwrite(data.frame(), self$out_ndt_all_file)
        data.table::fwrite(data.frame(), self$out_ndt_cvdonly_file)
        return(invisible(TRUE))
      }

      # standardize key columns
      ot_safety <- ot_safety %>%
        dplyr::rename(name_drug = drug_name, chemblId = chembl_id) %>%
        dplyr::select(-any_of("Gene"))

      # lightly parse indication field to a compact token (kept from original code)
      if ("indications" %in% names(ot_safety)) {
        ot_safety$indications <- sub(".*(efoName': '[^']+'?).*", "\\1", ot_safety$indications)
      } else {
        ot_safety$indications <- NA_character_
      }

      # ---------------------------- Merge / disease subset ----------------------------

      message("▶ Merging drug table with OT safety …")
      full_drug_safety <- safe_merge(
        drug_df,
        ot_safety,
        by = c("chemblId", "name_drug"),
        all.x = TRUE,
        allow.cartesian = TRUE
      )

      # Replace explicit “heart failure|cardiac failure” with the full disease_regex
      disease_term <- self$disease_regex
      full_drug_safety_hf_drugs <- dplyr::filter(
        full_drug_safety,
        grepl(disease_term, indications, ignore.case = TRUE)
      )

      message(sprintf("💾 Writing disease-like (by disease_regex) subset → %s", self$out_hf_file))
      data.table::fwrite(full_drug_safety_hf_drugs, self$out_hf_file)

      # ----------------------- Non-disease-term branch --------------------------

      nondiseaseterm_drug_df <- dplyr::filter(
        drug_df,
        Interacts_with_non_disease_terms_drug_only == "Yes"
      )

      nondiseaseterm_drug_safety <- safe_merge(
        nondiseaseterm_drug_df,
        ot_safety,
        by = c("chemblId", "name_drug"),
        all.x = TRUE,
        allow.cartesian = TRUE
      )

      # Use the *same* disease_regex to flag FDA adverse events that match disease terms
      if ("fda_adverse_drug_event" %in% names(nondiseaseterm_drug_safety)) {
        nondiseaseterm_drug_safety[, Cardiovascular_FDA_adverse_effect := ifelse(
          grepl(self$disease_regex, fda_adverse_drug_event, ignore.case = TRUE),
          "Yes", "No"
        )]
      } else {
        nondiseaseterm_drug_safety$Cardiovascular_FDA_adverse_effect <- NA_character_
        warning("⚠️ Column 'fda_adverse_drug_event' missing in safety table; cannot flag cardiovascular adverse effects.")
      }

      nondiseaseterm_drug_safety <- unique(nondiseaseterm_drug_safety)

      nondiseaseterm_drug_safety_cvd_side <- dplyr::filter(
        nondiseaseterm_drug_safety,
        Cardiovascular_FDA_adverse_effect == "Yes"
      )

      message(sprintf("💾 Writing nondiseaseterm all → %s", self$out_ndt_all_file))
      data.table::fwrite(nondiseaseterm_drug_safety, self$out_ndt_all_file)

      message(sprintf("💾 Writing nondiseaseterm CVD-only → %s", self$out_ndt_cvdonly_file))
      data.table::fwrite(nondiseaseterm_drug_safety_cvd_side, self$out_ndt_cvdonly_file)

      message("✅ Non-disease-term drug safety profiling complete.")
      invisible(TRUE)
    }
  )
)

# ----------------------------- Execute ------------------------------
runner <- NonDiseaseTermDrugProfile$new()
invisible(runner$run())
