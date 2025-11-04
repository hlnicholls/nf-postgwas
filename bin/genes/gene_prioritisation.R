#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(magrittr)
  library(here)
  library(data.table)
  library(R6)
})

# Shared pipeline config (defines: traits, var_file, locus_blocks_output, pops_output_path,
# coloc_output_path, enrichment_output_path, databases, prioritised_genes_path,
# var_path, hic_results_path (optional), custom_trait_name, disease_terms/regex, GTEx_tissue_terms/regex, etc.)
source("config_R.R")

# ----------------------------- Utilities -----------------------------

combine_values <- function(x) {
  ux <- unique(na.omit(x))
  if (length(ux) > 0) paste(ux, collapse = ",") else NA
}

clean_duplicates <- function(locus_name) {
  unique_genes <- unique(unlist(strsplit(locus_name, ",\\s*")))
  paste(unique_genes, collapse = ", ")
}

nzchar_or <- function(x) ifelse(is.na(x), "", x)  # avoid grepl on NAs

is_missing_or_empty_file <- function(path) {
  if (!file.exists(path)) return(TRUE)
  info <- tryCatch(file.info(path), error = function(e) NULL)
  if (is.null(info)) return(TRUE)
  is.na(info$size) || info$size == 0
}

# Helper: a data.table is "empty" if 0 rows OR all columns are length-0/NA
is_empty_dt <- function(dt) {
  is.null(dt) || nrow(dt) == 0 ||
    all(vapply(dt, function(col) length(col) == 0 || all(is.na(col)), logical(1)))
}

# Ensure columns exist; if missing, add with NA of the given mode
ensure_columns <- function(dt, cols, mode = "character") {
  for (cn in cols) {
    if (!(cn %in% names(dt))) {
      if (mode == "numeric") dt[[cn]] <- as.numeric(NA)
      else dt[[cn]] <- NA_character_
    }
  }
  dt
}

# ----------------------------- R6 Class ------------------------------

GenePrioritiser <- R6::R6Class(
  "GenePrioritiser",
  public = list(
    # ---- configuration (from config_R.R) ----
    traits = NULL,
    var_file = NULL,
    locus_blocks_output = NULL,
    pops_output_path = NULL,
    coloc_output_path = NULL,
    enrichment_output_path = NULL,
    databases = NULL,
    prioritised_genes_path = NULL,
    var_path = NULL,
    hic_results_path = NULL,
    custom_trait_name = NULL,

    # keyword lists/regex from config
    disease_terms = NULL,
    disease_regex = NULL,
    GTEx_tissue_terms = NULL,
    tissue_regex = NULL,

    # ---- dynamic/derived ----
    HiC_path = NULL,
    custom_pp4_all_file = NULL,
    custom_all_variants_pp4_file = NULL,
    has_custom_coloc_flag = NULL,        # dynamic output column name: has_<CUSTOM>_Coloc_H4_>0.8

    initialize = function() {
      # Bring in config values
      self$traits                 <- get("traits",                 envir = .GlobalEnv)
      self$var_file               <- get("var_file",               envir = .GlobalEnv)
      self$locus_blocks_output    <- get("locus_blocks_output",    envir = .GlobalEnv)
      self$pops_output_path       <- get("pops_output_path",       envir = .GlobalEnv)
      self$coloc_output_path      <- get("coloc_output_path",      envir = .GlobalEnv)
      self$enrichment_output_path <- get("enrichment_output_path", envir = .GlobalEnv)
      self$databases              <- get("databases",              envir = .GlobalEnv)
      self$prioritised_genes_path <- get("prioritised_genes_path", envir = .GlobalEnv)
      self$var_path               <- if (exists("var_path", envir = .GlobalEnv)) get("var_path", envir = .GlobalEnv) else NULL
      self$hic_results_path       <- if (exists("hic_results_path", envir = .GlobalEnv)) get("hic_results_path", envir = .GlobalEnv) else NULL
      self$custom_trait_name      <- get("custom_trait_name",      envir = .GlobalEnv)

      # disease/tissue terms & regex (compiled in config_R.R by the shim)
      self$disease_terms <- if (exists("disease_terms", envir = .GlobalEnv)) get("disease_terms", envir = .GlobalEnv) else character()
      self$GTEx_tissue_terms  <- if (exists("GTEx_tissue_terms",  envir = .GlobalEnv)) get("GTEx_tissue_terms",  envir = .GlobalEnv) else character()

      # Fallback regex defaults if shim didn't define compiled versions
      default_disease <- "cardio|cardiac|atrial|myocard|arrhyt|vascular|heart|hypertroph|dilated|left ventric|right ventric|obesity|obes"
      default_tissue  <- "heart|artery|aorta|atrial|adipose|ventricle|coronary"

      self$disease_regex <- if (exists("disease_regex", envir = .GlobalEnv)) get("disease_regex", envir = .GlobalEnv) else default_disease
      self$tissue_regex  <- if (exists("tissue_regex",  envir = .GlobalEnv)) get("tissue_regex",  envir = .GlobalEnv)  else default_tissue

      # Choose Hi-C file: prefer configured, else fallback under var_path
      self$HiC_path <- if (!is.null(self$hic_results_path)) {
        self$hic_results_path
      } else if (!is.null(self$var_path)) {
        file.path(self$var_path, "all_traits_in_ld_HiC.tsv")
      } else {
        "all_traits_in_ld_HiC.tsv"
      }

      # Custom coloc output files (produced by your custom coloc processor)
      self$custom_pp4_all_file          <- file.path(self$coloc_output_path, paste0(self$custom_trait_name, "_coloc_pp4_all.csv"))
      self$custom_all_variants_pp4_file <- file.path(self$coloc_output_path, paste0(self$custom_trait_name, "_coloc_all_variants_pp4.csv"))

      # Dynamic output column name (replaces `has_HF_Coloc_H4_>0.8`)
      self$has_custom_coloc_flag <- paste0("has_", self$custom_trait_name, "_Coloc_H4_>0.8")
    },

    run = function() {
      # 1) Base gene list from var_file + locus grouping
      genes <- self$build_gene_table_()

      # 2) Collapsed per gene
      collapsed_df <- self$collapse_by_gene_(genes)

      # 3) Merge PoPS flag (empty → all "No")
      collapsed_df <- self$attach_pops_(collapsed_df)

      # 4) Merge CUSTOM coloc (PP4≥0.8) per gene (empty → all "No")
      collapsed_df <- self$attach_custom_coloc_flag_(collapsed_df)

      # 5) EnrichR (most recent) — if absent, create NA columns and merge
      collapsed_df <- self$attach_enrichr_(collapsed_df)

      # 6) Optional Hi-C annotation — if absent, create NA HiC_tissue and merge
      collapsed_df <- self$attach_hic_(collapsed_df)

      # 7) GTEx eQTL coloc PP4≥0.8 merged across traits
      merged_df <- self$attach_gtex_coloc_(collapsed_df)

      # 8) Condense by gene
      condensed_df <- self$condense_per_gene_(merged_df)

      # 9) ClinGen (classification + disease)
      condensed_df <- self$attach_clingen_(condensed_df)

      # 10) DGIdb approved drugs list (with interaction score)
      condensed_df <- self$attach_dgidb_(condensed_df)

      # 11) IMPC phenotypes
      condensed_df <- self$attach_impc_(condensed_df)

      # 12) Compute prioritisation scores (uses custom coloc & keyword regexes)
      output <- self$compute_prioritisation_scores_(condensed_df)

      # 13) Curate fields to RV publication-like table + OpenTargets, DGIdb category
      dt_output <- self$curate_publication_table_(output)

      # 14) Final cleaning and write
      dt_output <- self$final_cleaning_(dt_output)
      dir.create(self$prioritised_genes_path, recursive = TRUE, showWarnings = FALSE)
      fwrite(dt_output, file.path(self$prioritised_genes_path, 'Prioritised_genes.csv'))

      cat("Gene prioritisation completed successfully\n")
      invisible(TRUE)
    },

    # ----------------------------- Steps ------------------------------

    build_gene_table_ = function() {
      genes <- data.table::fread(self$var_file) %>%
        dplyr::select(Gene_Symbol, Phenotype, Locus_name, lead_snp) %>%
        filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
        rename(Nearest_Gene_10kb = Gene_Symbol)

      locus_groups <- data.table::fread(self$locus_blocks_output, select = c('Lead_SNP', 'Locus_number')) %>%
        rename(lead_snp = Lead_SNP) %>%
        unique()

      genes <- merge(genes, locus_groups, by = 'lead_snp', all.x = TRUE, allow.cartesian = TRUE) %>%
        dplyr::select(-lead_snp) %>%
        unique()

      genes
    },

    collapse_by_gene_ = function(genes) {
      collapsed_df <- genes[, c(
        list(Phenotype = paste(unique(Phenotype), collapse = ", ")),
        lapply(.SD, function(x) paste(unique(x), collapse = ", "))
      ),
      by = Nearest_Gene_10kb,
      .SDcols = setdiff(names(genes), c("Nearest_Gene_10kb", "Phenotype"))]

      collapsed_df %>% filter(!is.na(Nearest_Gene_10kb))
    },

    attach_pops_ = function(df) {
      pops_path <- file.path(self$pops_output_path, 'gwas_all_loci_top_pops_genes.txt')

      if (is_missing_or_empty_file(pops_path)) {
        message(sprintf("PoPS file missing/empty: %s. Setting Top_PoPS_Score_per_Locus='No' for all genes.", pops_path))
        df$Top_PoPS_Score_per_Locus <- "No"
        return(df)
      }

      pops_data <- tryCatch(
        fread(pops_path, select = c('Nearest_Gene_10kb', 'Top_PoPS_Score_per_Locus')),
        error = function(e) data.table()
      )

      if (is_empty_dt(pops_data)) {
        message(sprintf("PoPS table is empty: %s. Setting Top_PoPS_Score_per_Locus='No' for all genes.", pops_path))
        df$Top_PoPS_Score_per_Locus <- "No"
        return(df)
      }

      df <- merge(df, pops_data, all.x = TRUE)
      df$Top_PoPS_Score_per_Locus[df$Top_PoPS_Score_per_Locus == "" | is.na(df$Top_PoPS_Score_per_Locus)] <- "No"
      df
    },

    attach_custom_coloc_flag_ = function(df) {
      # CUSTOM coloc PP4≥0.8 per gene
      if (is_missing_or_empty_file(self$custom_pp4_all_file)) {
        message(sprintf("Custom coloc PP4 file missing/empty: %s. Setting '%s' = 'No' for all genes.",
                        self$custom_pp4_all_file, self$has_custom_coloc_flag))
        df[[self$has_custom_coloc_flag]] <- "No"
        return(df)
      }

      coloc_data <- tryCatch(fread(self$custom_pp4_all_file), error = function(e) data.table())
      if (is_empty_dt(coloc_data)) {
        message(sprintf("Custom coloc PP4 table has 0 rows: %s. Setting '%s' = 'No' for all genes.",
                        self$custom_pp4_all_file, self$has_custom_coloc_flag))
        df[[self$has_custom_coloc_flag]] <- "No"
        return(df)
      }

      coloc_data <- coloc_data %>%
        dplyr::select(Nearest_Gene_10kb) %>%
        unique() %>%
        mutate(tmp_flag = "Yes")

      df <- merge(df, coloc_data, by = "Nearest_Gene_10kb", all.x = TRUE)
      df[[self$has_custom_coloc_flag]] <- ifelse(!is.na(df$tmp_flag) & df$tmp_flag == "Yes", "Yes", "No")
      df$tmp_flag <- NULL

      df %>% filter(!is.na(Nearest_Gene_10kb))
    },

    attach_enrichr_ = function(df) {
      # Expected EnrichR-derived columns used later in curate_publication_table_
      expected_cols <- c(
        "ClinVar_2019",
        "GO_Biological_Process_2023",
        "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023",
        "KEGG_2021_Human",
        "GTEx_Tissues_V8_2023",
        "GWAS_Catalog_2023",
        "OMIM_Expanded"
      )

      enrichr_files <- list.files(path = self$enrichment_output_path,
                                  pattern = "^All_genes_annotated_with_EnrichR_.*\\.csv$",
                                  full.names = TRUE)

      # Fallback: create NA columns for all genes and merge
      make_na_enrichr <- function(df) {
        na_tab <- data.table(Nearest_Gene_10kb = unique(df$Nearest_Gene_10kb))
        for (cn in expected_cols) na_tab[[cn]] <- NA_character_
        merge(df, na_tab, by = 'Nearest_Gene_10kb', all.x = TRUE)
      }

      if (length(enrichr_files) == 0) {
        message("No EnrichR results file found; creating NA EnrichR columns.")
        return(make_na_enrichr(df))
      }

      enrichr_file <- enrichr_files[which.max(file.info(enrichr_files)$mtime)]
      cat(sprintf("Loading EnrichR results from: %s\n", enrichr_file))

      if (is_missing_or_empty_file(enrichr_file)) {
        message("EnrichR file exists but is empty; creating NA EnrichR columns.")
        return(make_na_enrichr(df))
      }

      enrichr_results <- tryCatch(fread(enrichr_file), error = function(e) data.table())
      if (is_empty_dt(enrichr_results)) {
        message("EnrichR table is empty; creating NA EnrichR columns.")
        return(make_na_enrichr(df))
      }

      # First column is gene symbol; normalise its name
      enrichr_results <- enrichr_results %>%
        rename(Nearest_Gene_10kb = 1)

      # Ensure all expected columns exist (add as NA if missing)
      enrichr_results <- ensure_columns(enrichr_results, expected_cols, mode = "character")

      merge(df, enrichr_results, by = 'Nearest_Gene_10kb', all.x = TRUE)
    },

    attach_hic_ = function(df) {
      make_na_hic <- function(df) {
        na_tab <- data.table(Nearest_Gene_10kb = unique(df$Nearest_Gene_10kb),
                             HiC_tissue = NA_character_)
        merge(df, na_tab, by = 'Nearest_Gene_10kb', all.x = TRUE)
      }

      if (file.exists(self$HiC_path) && !is_missing_or_empty_file(self$HiC_path)) {
        hic_results <- tryCatch(fread(self$HiC_path), error = function(e) data.table())
        if (is_empty_dt(hic_results)) {
          message(sprintf("HiC file is empty at %s; creating NA HiC_tissue for all genes.", self$HiC_path))
          return(make_na_hic(df))
        }
        hic_results <- hic_results %>%
          dplyr::select(HiC_gene, HiC_tissue) %>%
          rename(Nearest_Gene_10kb = HiC_gene) %>%
          unique()

        hic_results <- hic_results[, c(
          list(HiC_tissue = paste(unique(HiC_tissue), collapse = ", ")),
          lapply(.SD, function(x) paste(unique(x), collapse = ", "))
        ),
        by = Nearest_Gene_10kb,
        .SDcols = setdiff(names(hic_results), c("Nearest_Gene_10kb", "HiC_tissue"))]

        merge(df, hic_results, by = 'Nearest_Gene_10kb', all.x = TRUE)
      } else {
        message(sprintf("HiC file not found or empty at %s; creating NA HiC_tissue for all genes.", self$HiC_path))
        make_na_hic(df)
      }
    },

    attach_gtex_coloc_ = function(df) {
      coloc_data_list <- list()

      for (trait in self$traits) {
        coloc_file <- file.path(self$coloc_output_path, paste0(trait, "_eQTL_COLOC.tsv"))
        if (file.exists(coloc_file) && !is_missing_or_empty_file(coloc_file)) {
          gtex_coloc <- tryCatch(fread(coloc_file), error = function(e) data.table())
          if (!is_empty_dt(gtex_coloc)) {
            gtex_coloc <- gtex_coloc %>%
              dplyr::select(Gene_Symbol, tissue, PP4) %>%
              rename(Nearest_Gene_10kb = Gene_Symbol) %>%
              filter(PP4 >= 80) %>%
              group_by(Nearest_Gene_10kb) %>%
              summarise(across(everything(), combine_values), .groups = "drop") %>%
              mutate(`has_GTEx_coloc_H4_>0.8` = "Yes")
            coloc_data_list[[trait]] <- gtex_coloc
          } else {
            message(sprintf("GTEx coloc file is empty for trait %s: %s", trait, coloc_file))
          }
        } else {
          message(paste("File for trait", trait, "not found or empty:", coloc_file))
        }
      }

      if (length(coloc_data_list) > 0) {
        all_gtex_coloc_08 <- dplyr::bind_rows(coloc_data_list)
        merge(df, all_gtex_coloc_08, by = 'Nearest_Gene_10kb', all.x = TRUE)
      } else df
    },

    condense_per_gene_ = function(df) {
      out <- df %>%
        group_by(Nearest_Gene_10kb) %>%
        summarise(across(everything(), combine_values), .groups = "drop") %>%
        filter(!is.na(Nearest_Gene_10kb))

      out
    },

    attach_clingen_ = function(df) {
      clingen_path <- file.path(self$databases, 'ClinGen/clingen_gene_disease_validities_31JAN2024.csv')
      if (is_missing_or_empty_file(clingen_path)) {
        message("ClinGen table missing/empty; continuing without ClinGen.")
        return(df)
      }
      clingen <- tryCatch(fread(clingen_path), error = function(e) data.table())
      if (is_empty_dt(clingen)) {
        message("ClinGen table is empty; continuing without ClinGen.")
        return(df)
      }

      clingen <- clingen %>%
        rename(
          Nearest_Gene_10kb = Gene,
          ClinGen_Disease = Disease,
          ClinGen_Classification = Classification
        ) %>%
        dplyr::select(Nearest_Gene_10kb, ClinGen_Disease, ClinGen_Classification) %>%
        unique() %>%
        group_by(Nearest_Gene_10kb) %>%
        summarise(across(everything(), combine_values), .groups = "drop")

      merge(df, clingen, by = 'Nearest_Gene_10kb', all.x = TRUE)
    },

    attach_dgidb_ = function(df) {
      dgidb_interact_path <- file.path(self$databases, 'dgidb/interactions.tsv')
      if (is_missing_or_empty_file(dgidb_interact_path)) {
        message("DGIdb interactions missing/empty; continuing without DGIdb interactions.")
        return(df)
      }
      dgidb_interact <- tryCatch(fread(dgidb_interact_path), error = function(e) data.table())
      if (is_empty_dt(dgidb_interact)) {
        message("DGIdb interactions are empty; continuing without DGIdb interactions.")
        return(df)
      }

      dgidb_interact <- dgidb_interact %>%
        rename(
          Nearest_Gene_10kb = gene_claim_name,
          dgidb_drug_name = drug_name,
          dgidb_approved_drug = approved,
          dgidb_interaction_score = interaction_score
        ) %>%
        dplyr::select(Nearest_Gene_10kb, dgidb_drug_name, dgidb_approved_drug, dgidb_interaction_score)

      dgidb_interact$dgidb_drug_name <- ifelse(
        is.na(dgidb_interact$dgidb_interaction_score),
        dgidb_interact$dgidb_drug_name,
        paste(dgidb_interact$dgidb_drug_name, sprintf("(%0.4f)", as.numeric(dgidb_interact$dgidb_interaction_score)))
      )
      dgidb_interact <- dgidb_interact[dgidb_interact$dgidb_drug_name != 'NULL (NA)', ] %>%
        filter(grepl("TRUE", dgidb_approved_drug)) %>%
        group_by(Nearest_Gene_10kb) %>%
        summarise(
          `DGIdb_drug(s) (interaction score)` = paste(dgidb_drug_name, collapse = ", "),
          dgidb_approved_drug = paste(dgidb_approved_drug, collapse = ", "),
          .groups = 'drop'
        ) %>%
        dplyr::select(Nearest_Gene_10kb, `DGIdb_drug(s) (interaction score)`)

      merge(df, dgidb_interact, by = 'Nearest_Gene_10kb', all.x = TRUE)
    },

    attach_impc_ = function(df) {
      impc_path <- file.path(self$databases, 'impc/phenotypeHitsPerGene_20250416.csv')
      if (is_missing_or_empty_file(impc_path)) {
        message("IMPC file missing/empty; continuing without IMPC annotations.")
        return(df)
      }
      impc <- tryCatch(fread(impc_path), error = function(e) data.table())
      if (is_empty_dt(impc)) {
        message("IMPC table is empty; continuing without IMPC annotations.")
        return(df)
      }

      impc <- impc %>%
        rename(
          Nearest_Gene_10kb = `Gene Symbol`,
          MGI_IMPC_Phenotypes = `Phenotype Hits`
        ) %>%
        mutate(Nearest_Gene_10kb = toupper(Nearest_Gene_10kb)) %>%
        dplyr::select(Nearest_Gene_10kb, MGI_IMPC_Phenotypes)

      df <- merge(df, impc, all.x = TRUE) %>% filter(!is.na(Nearest_Gene_10kb))
      df
    },

    compute_prioritisation_scores_ = function(condensed_df) {
      # Keyword-aware scoring using user-provided tissue & disease regexes
      output <- condensed_df %>%
        mutate(
          Pops = ifelse(Top_PoPS_Score_per_Locus == 'Yes', 1, 0),
          HiC_or_GTEx = ifelse(
            !is.na(HiC_tissue) |
              grepl(self$tissue_regex, nzchar_or(GTEx_Tissues_V8_2023), ignore.case = TRUE),
            1, 0
          ),
          MGI_or_OMIM_or_ClinGen =
            ifelse(
              grepl(self$disease_regex, nzchar_or(MGI_IMPC_Phenotypes), ignore.case = TRUE) |
                grepl(self$disease_regex, nzchar_or(ClinGen_Disease),      ignore.case = TRUE) |
                grepl(self$disease_regex, nzchar_or(OMIM_Expanded),        ignore.case = TRUE),
              1, 0
            )
        )

      # Start from CUSTOM PP4 gene flag (Yes/No) — if missing or empty earlier, it was set to "No".
      if (!(self$has_custom_coloc_flag %in% names(output))) {
        output[[self$has_custom_coloc_flag]] <- "No"
      }
      output[[self$has_custom_coloc_flag]] <- ifelse(grepl("Yes", nzchar_or(output[[self$has_custom_coloc_flag]])), "Yes", "No")
      output <- output %>%
        mutate(HF_coloc = ifelse(.data[[self$has_custom_coloc_flag]] == "Yes", 1, 0)) # helper for score sum

      cols_of_interest <- c('Pops', 'HiC_or_GTEx', 'MGI_or_OMIM_or_ClinGen', 'HF_coloc')
      for (c in cols_of_interest) {
        cat(stringr::str_interp("${c} ${sum(output[[c]], na.rm = TRUE)}\n"))
      }

      output <- output %>%
        mutate(Gene_Prioritisation_Score = rowSums(dplyr::select(., all_of(cols_of_interest)), na.rm = TRUE)) %>%
        mutate(across(where(is.character), ~ sub("^,", "", .)))

      # ---- Additional +1 boost using all-variants CUSTOM PP4 table
      if (is_missing_or_empty_file(self$custom_all_variants_pp4_file)) {
        warning(sprintf("Custom coloc all-variants PP4 file missing/empty: %s. Creating zero-boost flag for all genes.", self$custom_all_variants_pp4_file))
        output$has_HF_Coloc_H4_.0.8_new <- "No"
      } else {
        coloc_tbl <- tryCatch(fread(self$custom_all_variants_pp4_file), error = function(e) data.table())
        if (is_empty_dt(coloc_tbl)) {
          warning(sprintf("Custom coloc all-variants PP4 table has 0 rows: %s. Creating zero-boost flag for all genes.", self$custom_all_variants_pp4_file))
          output$has_HF_Coloc_H4_.0.8_new <- "No"
        } else {
          output_uncollapse <- output %>%
            dplyr::select(all_of(c("Phenotype", "Nearest_Gene_10kb", "Locus_name"))) %>%
            separate_longer_delim(Locus_name, delim = ', ') %>%
            mutate(Locus_number = suppressWarnings(as.integer(Locus_name)))

          output_uncollapse[[paste0("has_", self$custom_trait_name, "_Coloc_H4_.0.8_new")]] <-
            ifelse(output_uncollapse$Locus_name %in% coloc_tbl$Locus_name, 'Yes', 'No')

          new_custom_coloc_genes <- output_uncollapse %>%
            filter(.data[[paste0("has_", self$custom_trait_name, "_Coloc_H4_.0.8_new")]] == "Yes") %>%
            distinct(Nearest_Gene_10kb) %>%
            pull(Nearest_Gene_10kb)

          output <- output %>%
            mutate(has_HF_Coloc_H4_.0.8_new = ifelse(Nearest_Gene_10kb %in% new_custom_coloc_genes, 'Yes', 'No'))
        }
      }

      # Apply the zero/boost properly (+1 only where the new flag is Yes and previous HF_coloc==0)
      if ("has_HF_Coloc_H4_.0.8_new" %in% names(output)) {
        output <- output %>%
          mutate(Gene_Prioritisation_Score =
                   dplyr::case_when(
                     HF_coloc == 0 & has_HF_Coloc_H4_.0.8_new == 'Yes' ~ Gene_Prioritisation_Score + 1,
                     TRUE ~ Gene_Prioritisation_Score
                   ))
      }

      # Coalesce boosted flag into the dynamic flag if present (ensures a consistent single column downstream)
      if ("has_HF_Coloc_H4_.0.8_new" %in% names(output)) {
        target <- self$has_custom_coloc_flag
        newcol <- "has_HF_Coloc_H4_.0.8_new"
        if (target %in% names(output)) {
          output[[target]] <- dplyr::coalesce(output[[newcol]], output[[target]])
          output[[newcol]] <- NULL
        } else {
          output <- dplyr::rename(output, !!target := !!rlang::sym(newcol))
        }
      }

      # Drop the numeric helper
      output <- dplyr::select(output, -HF_coloc)
      output
    },

    curate_publication_table_ = function(output) {
      cols_core <- c(
        "Phenotype", "Nearest_Gene_10kb", "Locus_number", "Locus_name",
        "Gene_Prioritisation_Score",
        "DGIdb_drug(s) (interaction score)", "Top_PoPS_Score_per_Locus",
        "MGI_IMPC_Phenotypes", "ClinGen_Disease", 
        "OMIM_Expanded", "HiC_tissue", "ClinVar_2019",
        "GO_Biological_Process_2023", "GO_Cellular_Component_2023",
        "GO_Molecular_Function_2023", "KEGG_2021_Human",
        "GTEx_Tissues_V8_2023", "GWAS_Catalog_2023"
      )

      dt <- dplyr::select(output, all_of(cols_core))

      # OpenTargets merges
      ot_drugs_path <- file.path(self$databases, 'opentargets/OT_drug_interactions.csv')
      ot_warn_path  <- file.path(self$databases, 'opentargets/drugwarnings.csv')
      ot_pgx_path   <- file.path(self$databases, 'opentargets/pharmacogenomics.csv')

      if (file.exists(ot_drugs_path) && !is_missing_or_empty_file(ot_drugs_path) &&
          file.exists(ot_warn_path)  && !is_missing_or_empty_file(ot_warn_path)) {

        ot_drugs <- tryCatch(fread(ot_drugs_path), error = function(e) data.table())
        ot_warnings <- tryCatch(fread(ot_warn_path), error = function(e) data.table())

        if (!is_empty_dt(ot_drugs) && !is_empty_dt(ot_warnings)) {
          ot_warnings$chemblIds <- gsub("\\['|'\\]|'", "", ot_warnings$chemblIds)
          ot_warnings <- ot_warnings %>% tidyr::separate_rows(chemblIds, sep = ",\\s*") %>%
            dplyr::select(chemblIds, toxicityClass, country, description, efo_term)

          ot_dt <- merge(ot_drugs, ot_warnings, by = 'chemblIds', all.x = TRUE, allow.cartesian = TRUE) %>%
            rename(Nearest_Gene_10kb = name_mechanisms) %>%
            unique() %>%
            filter(!is.na(chemblIds) & chemblIds != "")

          if (file.exists(ot_pgx_path) && !is_missing_or_empty_file(ot_pgx_path)) {
            ot_pharmgkb <- tryCatch(fread(ot_pgx_path), error = function(e) data.table())
            if (!is_empty_dt(ot_pharmgkb)) {
              ot_pharmgkb <- ot_pharmgkb %>%
                rename(chemblIds = drugId) %>%
                dplyr::select(chemblIds, datasourceId, pgxCategory, phenotypeText, directionality, variantRsId, isDirectTarget)
              ot_dt <- merge(ot_dt, ot_pharmgkb, by = 'chemblIds', all.x = TRUE, allow.cartesian = TRUE)
            }
          }

          ot_dt <- ot_dt %>%
            group_by(Nearest_Gene_10kb) %>%
            summarise(across(everything(), combine_values), .groups = "drop") %>%
            dplyr::select(Nearest_Gene_10kb, name_drug)

          dt <- merge(dt, as.data.frame(ot_dt), by = 'Nearest_Gene_10kb', all.x = TRUE) %>%
            rename(`OT_drug(s)` = name_drug)
        } else {
          message("OpenTargets tables empty; skipping OT merges.")
        }
      } else {
        message("OpenTargets files missing/empty; skipping OT merges.")
      }

      dgidb_path <- file.path(self$databases, 'dgidb/categories.tsv')
      if (file.exists(dgidb_path) && !is_missing_or_empty_file(dgidb_path)) {
        dgidb <- tryCatch(fread(dgidb_path), error = function(e) data.table())
        if (!is_empty_dt(dgidb)) {
          dgidb <- dgidb %>%
            rename(Nearest_Gene_10kb = Gene) %>%
            mutate(Druggable_dgidb = ifelse(category == 'DRUGGABLE GENOME', 1, 0)) %>%
            dplyr::select(Nearest_Gene_10kb, Druggable_dgidb) %>%
            unique() %>%
            filter(Druggable_dgidb == 1)

          dt <- merge(dt, dgidb, by = 'Nearest_Gene_10kb', all.x = TRUE)
        } else {
          message("DGIdb categories empty; skipping druggable flag.")
        }
      } else {
        message("DGIdb categories missing/empty; skipping druggable flag.")
      }

      # Reorder and append the dynamic custom coloc column near the score
      dt <- dt %>%
        relocate(`OT_drug(s)`, .after = "GTEx_Tissues_V8_2023") %>%
        relocate(`DGIdb_drug(s) (interaction score)`, .after = `OT_drug(s)`) %>%
        relocate(Druggable_dgidb, .after = `DGIdb_drug(s) (interaction score)`)

      # Ensure presence of dynamic custom coloc flag; keep as Yes/No
      add_col <- self$has_custom_coloc_flag
      if (!(add_col %in% names(output))) {
        output[[add_col]] <- "No"
      }
      dt <- merge(dt,
                  output[, c("Nearest_Gene_10kb", add_col), drop = FALSE],
                  by = "Nearest_Gene_10kb", all.x = TRUE) %>%
        relocate(all_of(add_col), .after = "Gene_Prioritisation_Score")

      dt
    },

    final_cleaning_ = function(dt_output) {
      # Remove leading or trailing ", NA" artifacts across annotation cols
      if (ncol(dt_output) >= 8) {
        dt_output[, 8:ncol(dt_output)] <- lapply(dt_output[, 8:ncol(dt_output)], function(x) gsub("NA,\\s*", "", x))
        dt_output[, 8:ncol(dt_output)] <- lapply(dt_output[, 8:ncol(dt_output)], function(x) gsub(",\\s*NA", "", x))
      }

      dt_output <- dt_output %>% filter(Nearest_Gene_10kb != "")
      cat('Final check for duplicates (should say FALSE):\n')
      print(any(duplicated(dt_output$Nearest_Gene_10kb)))

      # Sort by locus number as before
      dt_output <- dt_output[order(as.numeric(dt_output$Locus_number)), ]
      dt_output
    }
  )
)

# ----------------------------- Execute ------------------------------

gp <- GenePrioritiser$new()
invisible(gp$run())
