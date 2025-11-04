#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(R6)
  library(data.table)
  library(dplyr)
})

#-----------------------------#
# Argument parsing
#-----------------------------#
opt <- OptionParser() |>
  add_option("--in",  type = "character", dest = "inp",  help = "Input <trait>_38_37.txt") |>
  add_option("--out", type = "character", dest = "out",  help = "Output <trait>_38_37_rsids.txt") |>
  add_option("--bim_rds", type = "character", dest = "bim_rds",
             help = "Path to prebuilt BIM RDS with columns: SNP, rsid_1kg") |>
  parse_args()

if (is.null(opt$inp) || is.null(opt$out) || is.null(opt$bim_rds)) {
  stop("Required: --in, --out, --bim_rds", call. = FALSE)
}

if (!grepl("^[A-Za-z0-9_]+_38_37\\.txt$", basename(opt$inp))) {
  stop("Input file name does not match expected <trait>_38_37.txt pattern: ", basename(opt$inp), call. = FALSE)
}

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
# Try to use all threads for data.table (same as original)
try(data.table::setDTthreads(percent = 100), silent = TRUE)

#-----------------------------#
# RsidMapper class (R6)
#-----------------------------#
RsidMapper <- R6Class(
  "RsidMapper",
  public = list(
    # Inputs / outputs
    inp = NULL,
    out = NULL,
    bim_rds = NULL,

    # Internal state
    bim_map = NULL,   # data.table with SNP, rsid_1kg
    snps = NULL,      # data.table read from input
    merged_data = NULL,

    initialize = function(inp, out, bim_rds) {
      self$inp <- inp
      self$out <- out
      self$bim_rds <- bim_rds
    },

    #' Run the full pipeline (same behavior as original)
    run = function() {
      self$load_bim_map()
      self$read_snps()
      self$ensure_snp_column()
      self$merge_with_bim()
      self$drop_plotting_qc_extras()
      self$write_output()
      invisible(self)
    },

    #' Load BIM RDS and keep only SNP, rsid_1kg (validate required columns)
    load_bim_map = function() {
      bim <- readRDS(self$bim_rds)
      req <- c("SNP", "rsid_1kg")
      miss <- setdiff(req, names(bim))
      if (length(miss)) {
        stop("BIM RDS missing columns: ", paste(miss, collapse = ", "))
      }
      self$bim_map <- as.data.table(bim)[, .(SNP, rsid_1kg)]
    },

    #' Read full input (preserve all columns, no type coercion beyond fread defaults)
    read_snps = function() {
      self$snps <- fread(self$inp, header = TRUE, showProgress = FALSE)
    },

    #' Ensure an SNP column exists, otherwise compose it as CHROM:GENPOS:ALLELE0:ALLELE1
    ensure_snp_column = function() {
      if (!("SNP" %in% names(self$snps))) {
        need <- c("CHROM", "GENPOS", "ALLELE0", "ALLELE1")
        miss <- setdiff(need, names(self$snps))
        if (length(miss)) {
          stop("Input missing SNP and required cols to construct it: ", paste(miss, collapse = ", "))
        }
        # NOTE: matches original behavior; NA values become "NA" string in paste
        self$snps[, SNP := paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep = ":")]
      }
    },

    #' Left-merge GWAS rows with BIM rsid map; keep ALL original columns + add rsid_1kg
    merge_with_bim = function() {
      self$merged_data <- merge(self$snps, self$bim_map, by = "SNP", all.x = TRUE)
    },

    #' Drop only plotting/QC extras if present (do NOT drop P/BETA/SE)
    drop_plotting_qc_extras = function() {
      drop_cols <- intersect(names(self$merged_data), c("INFO", "CHISQ", "LOG10P", "TEST"))
      if (length(drop_cols)) {
        self$merged_data <- dplyr::select(self$merged_data, -dplyr::all_of(drop_cols))
      }
    },

    #' Write final merged table as TSV (tab-separated, no quotes)
    write_output = function() {
      fwrite(self$merged_data, file = self$out, sep = "\t", quote = FALSE, showProgress = FALSE)
      message("Wrote: ", self$out)
    }
  )
)

# ----------------------------- Execute -----------------------------
RsidMapper$new(inp = opt$inp, out = opt$out, bim_rds = opt$bim_rds)$run()
