
#!/usr/bin/env Rscript
library(R6)
library(data.table)
library(here)

source("config_R.R")

LDSCLogCleaner <- R6Class(
  "LDSCLogCleaner",
  public = list(
    traits = NULL,
    ldsc_results_path = NULL,
    ldsc_log_path_for_table = NULL,

    initialize = function(traits, ldsc_results_path, ldsc_log_path_for_table) {
      self$traits <- traits
      self$ldsc_results_path <- ldsc_results_path
      self$ldsc_log_path_for_table <- ldsc_log_path_for_table
    },

    clean_trait_log = function(trait) {
      message(paste('Processing trait:', trait))
      log_file <- file.path(self$ldsc_results_path, paste0('genetic_correlation_ldsc_res_ukb_', trait, '.log'))
      log_contents <- readLines(log_file)

      start_marker <- "p1"
      end_marker <- "Analysis finished"
      table_start <- which(grepl(start_marker, log_contents))
      table_end <- which(grepl(end_marker, log_contents)) - 1
      table_contents <- log_contents[table_start:table_end]

      phenotypes <- self$traits
      path_to_remove <- self$ldsc_log_path_for_table
      cleaned_table <- sapply(table_contents, function(line) {
        for (phenotype in phenotypes) {
          line <- gsub(paste0(path_to_remove, phenotype, "_GWAS_37_corr.txt"), phenotype, line)
        }
        return(line)
      }, USE.NAMES = FALSE)

      output_file_path <- file.path(self$ldsc_results_path, paste0(trait, "_corr_table.txt"))
      table_string <- paste(cleaned_table, collapse = "\n")
      table_df <- read.table(text = table_string, header = TRUE, sep = "", fill = TRUE, strip.white = TRUE)
      write.table(table_df, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
      message(paste('Cleaned table saved for', trait, 'at', output_file_path))
    },

    clean_all_traits = function() {
      message('Cleaning LDSC log files...')
      for (trait in self$traits) {
        self$clean_trait_log(trait)
      }
    },

    combine_tables = function() {
      message('Combining all cleaned tables into a single file...')
      path <- self$ldsc_results_path
      file_list <- list.files(path, pattern = "_corr_table\\.txt$", full.names = TRUE)
      all_data <- data.table()
      for (file in file_list) {
        message(paste('Reading file:', file))
        dt <- fread(file)
        dt[, file_name := basename(file)]
        all_data <- rbind(all_data, dt, fill = TRUE)
      }
      if ("p1" %in% colnames(all_data)) {
        all_data[, p1 := gsub(".*/([^/]+)_GWAS_.*", "\\1", p1)]
      }
      if ("p2" %in% colnames(all_data)) {
        all_data[, p2 := gsub(".*/([^/]+)_GWAS_.*", "\\1", p2)]
      }
      cleaned_output_file <- file.path(self$ldsc_results_path, "all_pheno_corr_table.txt")
      fwrite(all_data, file = cleaned_output_file, sep = "\t", quote = FALSE)
      message(paste('Combined cleaned data saved at', cleaned_output_file))
    },

    run = function() {
      self$clean_all_traits()
      self$combine_tables()
    }
  )
)

# ----------------------------- Execute -----------------------------
cleaner <- LDSCLogCleaner$new(traits, ldsc_results_path, ldsc_log_path_for_table)
cleaner$run()
