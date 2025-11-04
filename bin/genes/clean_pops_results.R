#!/usr/bin/env Rscript

# Clean and annotate PoPS results for all traits

suppressPackageStartupMessages({
  library(data.table)
  library(biomaRt)
  library(here)
})

source("config_R.R")

fallback_file <- paste0(databases,"/pops/example/data/utils/gene_annot_jun10.txt")

# Function to connect to Ensembl
try_ensembl <- function(dataset = "hsapiens_gene_ensembl") {
  mirrors <- c("www.ensembl.org", "useast.ensembl.org", "uswest.ensembl.org", "asia.ensembl.org")
  for (mirror in mirrors) {
    message(sprintf("Trying Ensembl mirror: %s", mirror))
    ensembl <- tryCatch({
      useMart("ensembl", dataset = dataset, host = paste0("https://", mirror))
    }, error = function(e) {
      message(sprintf("Failed to connect to Ensembl mirror: %s.", mirror))
      return(NULL)
    })
    if (!is.null(ensembl)) {
      message(sprintf("Connected to Ensembl mirror: %s", mirror))
      return(ensembl)
    }
  }
  message("All Ensembl mirrors are currently unavailable.")
  return(NULL)  # Return NULL if all mirrors fail
}

# Function to get gene information
get_gene_info <- function(ensembl_ids, ensembl) {
  if (!is.null(ensembl)) {
    gene_info <- tryCatch({
      getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
            filters = 'ensembl_gene_id',
            values = ensembl_ids,
            mart = ensembl)
    }, error = function(e) {
      message("Error querying Ensembl. Switching to local annotation file.")
      return(NULL)
    })
    if (!is.null(gene_info) && nrow(gene_info) > 0) {
      return(gene_info)
    }
  }

  # Fallback to local file
  message("Using local gene annotation file for mapping.")
  annotation <- fread(fallback_file)
  if (!all(c("ENSGID", "NAME") %in% colnames(annotation))) {
    stop("Fallback file does not have required columns: ENSGID and NAME.")
  }
  gene_info <- annotation[, .(ensembl_gene_id = ENSGID, hgnc_symbol = NAME)]
  return(gene_info)
}

# Main script logic
process_pops_files <- function(traits, pops_output_path, output_file_name) {
  pops_list <- list()
  
  for (trait in traits) {
    file_name <- paste0(trait, "_pops_all_features.preds")
    full_path <- file.path(pops_output_path, file_name)
    dt <- fread(full_path, sep = "\t", header = TRUE)
    setnames(dt, "PoPS_Score", paste0(trait, "_pops"))
    colnames(dt)[colnames(dt) != "ENSGID"] <- paste0(colnames(dt)[colnames(dt) != "ENSGID"], "_", trait)
    pops_list[[trait]] <- dt
  }
  
  output <- Reduce(function(x, y) merge(x, y, by = "ENSGID", all = TRUE), pops_list)
  selected_columns <- c("ENSGID", grep("pops", names(output), value = TRUE))
  output_selected <- output[, ..selected_columns]
  
  # Attempt to connect to Ensembl
  ensembl <- try_ensembl()
  
  ensembl_ids <- unique(output_selected$ENSGID)  # Extract unique ENSG IDs
  gene_info <- get_gene_info(ensembl_ids, ensembl)
  
  output_annotated <- merge(output_selected, gene_info, by.x = "ENSGID", by.y = "ensembl_gene_id", all.x = TRUE)
  setcolorder(output_annotated, c("hgnc_symbol", setdiff(names(output_annotated), "hgnc_symbol")))
  
  fwrite(output_annotated, paste0(pops_output_path, '/', output_file_name))
  message(sprintf("File saved: %s", paste0(pops_output_path, '/', output_file_name)))
}

# Process traits
process_pops_files(traits, pops_output_path, "pops_results_all_features_cleaned.csv")

cat("PoPS results cleaning completed successfully\n")
