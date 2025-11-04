#!/usr/bin/env Rscript

# Wakefield Approximate Bayes Factor fine-mapping
# - Uses canonical BF_1:0 formula from Wakefield
# - Builds 95% and 99% credible sets per (Phenotype, Locus_n) group
# - Output: finemapping_path/All_traits_credible_sets.txt (tab-separated)

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table)
  library(corrcoverage)
  library(R6)
})


source("config_R.R")

# ----------------------------------------------------------------------------------------------
# Wakefield fine-mapping
# ----------------------------------------------------------------------------------------------
WakefieldFineMapper <- R6::R6Class(
  classname = "WakefieldFineMapper",
  public = list(
    # Public fields
    wake = NULL,                # prior variance (W)
    loci = NULL,                # raw input loci (data.table)
    traits = NULL,              # vector of Phenotype names
    sample_sizes = NULL,        # numeric vector aligned with traits
    finemapping_path = NULL,    # output directory
    group_keys = NULL,          # keys used for grouping loci

    initialize = function(wake = 0.04,
                          traits,
                          sample_sizes,
                          finemapping_path) {
      self$wake             <- wake
      self$traits           <- traits
      self$sample_sizes     <- sample_sizes
      self$finemapping_path <- finemapping_path
      dir.create(self$finemapping_path, showWarnings = FALSE, recursive = TRUE)
    },

    load_loci = function(var_file) {
      stopifnot(file.exists(var_file))
      loci_dt <- data.table::fread(var_file, sep = "\t", header = TRUE, na.strings = c("", "NA"))
      if (!nrow(loci_dt)) stop("Input loci file has 0 rows: ", var_file)

      # Remove SNPs that are in LD but not associated with any traits (keeps CHROM present)
      loci_dt <- loci_dt[!is.na(CHROM)]

      # Distinct by (Phenotype, SNP) at minimum, but keep other columns
      loci_dt <- loci_dt %>% dplyr::distinct(Phenotype, SNP, .keep_all = TRUE)

      # Decide grouping keys: prefer Locus_n; fallback to lead_snp; fallback to Locus_name
      keys <- c()
      if ("Locus_n" %in% names(loci_dt)) {
        keys <- c("Phenotype", "Locus_n")
      } else if ("lead_snp" %in% names(loci_dt)) {
        keys <- c("Phenotype", "lead_snp")
        message("Grouping by (Phenotype, lead_snp) because Locus_n is absent.")
      } else if ("Locus_name" %in% names(loci_dt)) {
        keys <- c("Phenotype", "Locus_name")
        message("Grouping by (Phenotype, Locus_name) because Locus_n/lead_snp are absent.")
      } else {
        stop("Cannot determine grouping key: expected one of Locus_n, lead_snp, or Locus_name.")
      }

      self$group_keys <- keys
      self$loci <- loci_dt
      invisible(self)
    },

    add_sample_sizes = function() {
      if (is.null(self$loci)) stop("Call load_loci() first.")
      # Build a named map Phenotype -> N
      if (length(self$traits) != length(self$sample_sizes))
        stop("traits and sample_sizes lengths differ.")

      sample_size_map <- stats::setNames(self$sample_sizes, self$traits)

      # Attach N and compute Z, V, ABF using canonical BF_1:0
      self$loci <- self$loci %>%
        dplyr::mutate(
          N  = sample_size_map[Phenotype],
          Z  = BETA / SE,           # sign preserved
          V  = SE * SE
        ) %>%
        dplyr::mutate(
          ABF = sqrt(V / (V + self$wake)) * exp(0.5 * (Z * Z) * (self$wake / (V + self$wake)))
        )

      invisible(self)
    },

    # Compute PP within each (Phenotype, Locus) group
    compute_pp = function() {
      if (is.null(self$loci)) stop("Call load_loci() first.")
      if (is.null(self$group_keys)) stop("Grouping keys are not set.")

      # Split into groups and normalize ABF to PP per group
      group_syms <- rlang::syms(self$group_keys)
      split_groups <- self$loci %>%
        dplyr::group_by(!!group_syms[[1]], !!!group_syms[-1]) %>%
        dplyr::group_split(.keep = TRUE)

      # Normalize within group
      split_groups <- purrr::map(split_groups, ~ .x %>% dplyr::mutate(PP = ABF / sum(ABF, na.rm = TRUE)))

      # Recombine
      self$loci <- data.table::as.data.table(dplyr::bind_rows(split_groups))
      invisible(self)
    },

    # Build credible sets at the requested threshold (e.g., 0.95 or 0.99)
    # Returns a data.table with the *subset* of rows that are in the credible set for each group,
    # with a flag column name supplied (e.g., "credible_set" or "credible_set_95") set to "Yes".
    build_credible_set = function(threshold = 0.99, flag_col = "credible_set") {
      if (is.null(self$loci)) stop("Call load_loci() first.")
      if (is.null(self$group_keys)) stop("Grouping keys are not set.")
      if (!"PP" %in% names(self$loci)) stop("Call compute_pp() before building credible sets.")

      group_syms <- rlang::syms(self$group_keys)
      groups <- self$loci %>%
        dplyr::group_by(!!group_syms[[1]], !!!group_syms[-1]) %>%
        dplyr::group_split(.keep = TRUE)

      cs_rows <- purrr::map(groups, function(g) {
        # Order by decreasing PP
        g2 <- g %>% dplyr::arrange(dplyr::desc(PP)) %>% dplyr::mutate(cum_PP = cumsum(PP))

        cs <- tryCatch(
          corrcoverage::credset(g2$PP, thr = threshold),
          error = function(e) NULL
        )

        if (is.null(cs) || is.na(cs$credset)) {
          # No credible set; return empty (we'll mark "No" later via join)
          return(g2[0, ])
        }

        # Keep only rows in the credible set and mark
        out <- g2 %>% dplyr::slice(cs$credset)
        out[[flag_col]] <- "Yes"
        out
      })

      data.table::as.data.table(dplyr::bind_rows(cs_rows))
    },

    # Run the entire pipeline and write output
    run = function() {
      if (is.null(self$loci)) stop("Call load_loci() first.")
      message("Wakefield fine-mapping on groups: ", paste(self$group_keys, collapse = ", "))

      self$add_sample_sizes()
      self$compute_pp()

      # Build 99% and 95% credible sets (subsets with flags)
      cs99 <- self$build_credible_set(threshold = 0.99, flag_col = "credible_set")
      cs95 <- self$build_credible_set(threshold = 0.95, flag_col = "credible_set_95")

      # Merge flags back to full loci (outer join by keys + SNP identity)
      key_cols <- unique(c(self$group_keys, "SNP", "Phenotype", "CHROM", "lead_snp"))
      key_cols <- key_cols[key_cols %in% names(self$loci)]  # keep only existing

      loci_out <- self$loci %>%
        dplyr::left_join(cs99 %>% dplyr::select(dplyr::any_of(c(key_cols, "credible_set"))),
                         by = intersect(key_cols, names(cs99))) %>%
        dplyr::left_join(cs95 %>% dplyr::select(dplyr::any_of(c(key_cols, "credible_set_95"))),
                         by = intersect(key_cols, names(cs95))) %>%
        dplyr::mutate(
          credible_set     = ifelse(is.na(.data$credible_set), "No",  "Yes"),
          credible_set_95  = ifelse(is.na(.data$credible_set_95), "No", "Yes")
        )

      # Ensure uniqueness and original column order + new columns near front
      front_cols <- unique(c("Phenotype", self$group_keys, "credible_set", "credible_set_95"))
      front_cols <- front_cols[front_cols %in% names(loci_out)]
      all_cols   <- unique(c(front_cols, names(loci_out)))

      # ---- FIX: data.table column selection must use .. for character vectors
      loci_out_dt <- data.table::as.data.table(loci_out)
      loci_final  <- unique(loci_out_dt[, ..all_cols])

      # Write output (same location/name as your previous script)
      out_file <- file.path(self$finemapping_path, "All_traits_credible_sets.txt")
      data.table::fwrite(loci_final, out_file, sep = "\t")

      message("Wrote: ", out_file)
      invisible(loci_final)
    }
  )
)

# ----------------------------- Execute -----------------------------

wake_prior <- 0.04  # same default as your previous script

mapper <- WakefieldFineMapper$new(
  wake             = wake_prior,
  traits           = traits,
  sample_sizes     = sample_sizes,
  finemapping_path = finemapping_path
)

# Load input, run, and write outputs
mapper$load_loci(var_file = var_file)
mapper$run()
