#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(R6)
})

# ------------------------------------------------------------
# Load pipeline config (defines: output_path, databases,
# prioritised_genes, druggability_output_path, druggability_results,
# and disease_terms / disease_regex injected by the shim)
# ------------------------------------------------------------
source("config_R.R")

# ----------------------------- Utilities -----------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b
nzchar_or <- function(x) ifelse(is.na(x), "", x)

combine_values <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) > 0) paste(ux, collapse = ", ") else NA_character_
}

# -------------------------- R6 Drug Table ----------------------------

DrugTable <- R6::R6Class(
  "DrugTable",
  public = list(
    # Config from config_R.R
    output_path = NULL,
    databases = NULL,
    prioritised_genes = NULL,
    druggability_output_path = NULL,
    druggability_results = NULL,  # alias for backward compat
    disease_terms = NULL,
    disease_regex = NULL,

    # Data slots
    impc = NULL,
    priors = NULL,
    ot_drugs = NULL,
    ot_warnings = NULL,
    ot_pgx = NULL,

    initialize = function() {
      # Bind config values
      self$output_path            <- get("output_path",            envir = .GlobalEnv)
      self$databases              <- get("databases",              envir = .GlobalEnv)
      self$prioritised_genes      <- get("prioritised_genes",      envir = .GlobalEnv)
      self$druggability_output_path <- get("druggability_output_path", envir = .GlobalEnv)
      # Allow either alias to exist
      self$druggability_results   <- if (exists("druggability_results", envir = .GlobalEnv)) {
        get("druggability_results", envir = .GlobalEnv)
      } else {
        self$druggability_output_path
      }

      # disease terms/regex from shim (fallback to a sane default if not present)
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

      # Ensure output dir exists
      dir.create(self$druggability_output_path, showWarnings = FALSE, recursive = TRUE)
    },

    run = function() {
      self$load_inputs_()
      merged <- self$merge_opentargets_()
      prior  <- self$restrict_to_prioritised_(merged)
      self$emit_counts_(prior)

      # Determine disease / non-disease flags using user-provided regex
      disease_genes <- self$genes_with_disease_terms_(prior)
      self$write_disease_gene_list_(disease_genes)

      non_disease_genes <- self$genes_without_disease_terms_(prior)
      final_tbl <- self$final_output_(merged, non_disease_genes)

      # Write the final druggability table
      fwrite(final_tbl, file.path(self$druggability_output_path, "Druggability_table.csv"))

      invisible(TRUE)
    },

    load_inputs_ = function() {
      # prioritised genes (required)
      self$priors <- fread(self$prioritised_genes) %>%
        dplyr::select(Nearest_Gene_10kb, Phenotype, Gene_Prioritisation_Score) %>%
        mutate(Reported_Locus = NA_character_) %>%     # ensure column exists (back-compat)
        rename(Gene = Nearest_Gene_10kb) %>%
        filter(Gene != "")

      # OpenTargets files
      self$ot_drugs <- fread(file.path(self$databases, "opentargets/OT_drug_interactions.csv")) %>%
        rename(Gene = name_mechanisms)

      self$ot_warnings <- fread(file.path(self$databases, "opentargets/drugwarnings.csv")) %>%
        mutate(chemblIds = gsub("\\['|'\\]|'", "", chemblIds)) %>%
        separate_rows(chemblIds, sep = ",\\s*") %>%
        dplyr::select(chemblIds, toxicityClass, country, description, efo_term)

      self$ot_pgx <- fread(file.path(self$databases, "opentargets/pharmacogenomics.csv")) %>%
        rename(chemblIds = drugId) %>%
        dplyr::select(chemblIds, datasourceId, pgxCategory, phenotypeText, directionality, variantRsId, isDirectTarget)
    },

    merge_opentargets_ = function() {
      # priors + drugs
      df <- merge(self$priors, self$ot_drugs, by = "Gene", all.x = TRUE)

      # + warnings
      dt <- merge(df, self$ot_warnings, by = "chemblIds", all.x = TRUE)

      # Column order (retain the existing columns used later)
      dt <- dt %>%
        dplyr::select(
          Phenotype, Gene, Reported_Locus, Gene_Prioritisation_Score, chemblIds,
          toxicityClass, country, description, efo_term,
          name_drug, actionType, mechanismOfAction, targetName, terms_drug,
          targetType, description_mechanisms, entity_mechanisms, category_mechanisms,
          description_drug, category_drug, indications, approvedIndications
        ) %>%
        unique() %>%
        filter(!is.na(chemblIds) & chemblIds != "")

      # + PharmGKB
      d <- merge(dt, self$ot_pgx, by = "chemblIds", all.x = TRUE) %>%
        unique()

      d
    },

    restrict_to_prioritised_ = function(merged_dt) {
      # Keep only rows where the gene had a non-zero prioritisation score
      merged_dt %>% filter(Gene_Prioritisation_Score != 0)
    },

    emit_counts_ = function(prior_dt) {
      # total distinct genes with any drug interactions
      gene_drugs_count <- unique(na.omit(prior_dt$Gene))
      cat("[1] \"Number of prioritsed gene-drug interactions\"\n")
      cat(sprintf("[1] %d\n", length(gene_drugs_count)))
    },

    genes_with_disease_terms_ = function(prior_dt) {
      # Determine which prioritised genes have any drug whose terms match disease_regex
      if (!("terms_drug" %in% names(prior_dt))) {
        # If the upstream file layout is unexpected, treat as no matches
        return(character(0))
      }

      disease_regex <- self$disease_regex %||% ""
      if (isTRUE(is.na(disease_regex)) || nchar(disease_regex) == 0) {
        return(character(0))
      }

      disease_gene_drugs <- prior_dt %>%
        group_by(Gene) %>%
        mutate(
          .flag = any(grepl(disease_regex, tolower(nzchar_or(terms_drug)), ignore.case = TRUE))
        ) %>%
        filter(.flag) %>%
        dplyr::select(-.flag) %>%
        summarise(across(everything(), combine_values), .groups = "drop")

      unique(disease_gene_drugs$Gene)
    },

    write_disease_gene_list_ = function(genes_vec) {
      # Write disease-matched gene list as before
      out_df <- data.frame(Gene = genes_vec) %>%
        mutate(Disease_Drug_Interaction = "Yes")

      fwrite(out_df, file.path(self$druggability_output_path, "disease_drug_genes.txt"))

      cat("[1] \"Number of prioritsed gene-drug interactions\"\n")
      cat(sprintf("[1] %d\n", nrow(out_df)))
    },

    genes_without_disease_terms_ = function(prior_dt) {
      # Genes that only map to non-disease-term drugs (used to annotate the final table)
      if (!("terms_drug" %in% names(prior_dt))) return(character(0))

      disease_regex <- self$disease_regex %||% ""
      if (isTRUE(is.na(disease_regex)) || nchar(disease_regex) == 0) {
        # If no regex provided, no gene can be flagged as non-disease-only
        return(character(0))
      }

      nondisease_gene_drugs <- prior_dt %>%
        group_by(Gene) %>%
        mutate(.flag = any(grepl(disease_regex, tolower(nzchar_or(terms_drug)), ignore.case = TRUE))) %>%
        # keep those that DO NOT have any disease-term drug
        filter(!.flag) %>%
        dplyr::select(-.flag) %>%
        summarise(across(everything(), combine_values), .groups = "drop")

      unique(nondisease_gene_drugs$Gene)
    },

    final_output_ = function(merged_dt, non_disease_genes) {
      # Add Interacts_with_non_disease_terms_drug_only flag
      merged_dt$Interacts_with_non_disease_terms_drug_only <-
        ifelse(merged_dt$Gene %in% non_disease_genes, "Yes", "No")

      out <- merged_dt %>%
        dplyr::select(
          Phenotype, Gene, Reported_Locus, Gene_Prioritisation_Score,
          Interacts_with_non_disease_terms_drug_only,
          toxicityClass, country, chemblIds, name_drug,
          pgxCategory, phenotypeText, directionality, variantRsId, isDirectTarget
        )

      data.table::setnames(
        out,
        old = c("chemblIds", "pgxCategory", "phenotypeText", "directionality", "variantRsId", "isDirectTarget"),
        new = c("chemblId", "pharmGKB_pgxCategory", "pharmGKB_phenotype", "pharmGKB_directionality", "pharmGKB_variantRsId", "pharmGKB_isDirectTarget")
      )

      out
    }
  )
)

# ----------------------------- Execute ------------------------------

dt <- DrugTable$new()
invisible(dt$run())
