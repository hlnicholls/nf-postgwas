#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(R6)
  library(data.table)
  library(Sushi)
  library(fields)
  library(biomaRt)
  library(here)
  library(dplyr)
})

source("config_R.R")

Sys.setenv(CURL_CA_BUNDLE = "/etc/ssl/certs/ca-certificates.crt")

LocusZoomPlotter <- R6::R6Class(
  "LocusZoomPlotter",
  public = list(
    # config
    gwas_dir = NULL,
    plot_dir = NULL,
    leadVar_path = NULL,
    ld_path = NULL,
    genetic_map_path = NULL,
    phenotypes = NULL,
    plot_win = 500e3,   # +/- 500kb (1Mb window total)
    ensembl = NULL,

    # constructor: set config, resolve helper script, init ensembl, ensure dirs
    initialize = function(gwas_dir, plot_dir, leadVar_path, ld_path, genetic_map_path, phenotypes, plot_win = 500e3) {
      self$gwas_dir         <- gwas_dir
      self$plot_dir         <- plot_dir
      self$leadVar_path     <- leadVar_path
      self$ld_path          <- ld_path
      self$genetic_map_path <- genetic_map_path
      self$phenotypes       <- phenotypes
      self$plot_win         <- plot_win

      private$source_makeRegionalPlot()

      # Be tolerant to transient biomart hiccups; allow plots to proceed even if NULL
      self$ensembl <- tryCatch(
        useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl"),
        error = function(e) {
          message("Warning: biomaRt connection failed. Proceeding without gene track. Error: ", conditionMessage(e))
          NULL
        }
      )

      if (!dir.exists(self$plot_dir)) dir.create(self$plot_dir, recursive = TRUE)
    },

    run = function() {
      # mimic original: setwd(plot_dir)
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(self$plot_dir)

      for (pheno in self$phenotypes) {
        message(sprintf("Print regional plots for %s", pheno))

        # Wrap the entire phenotype block so one bad phenotype doesn't kill the run
        tryCatch({
          message("------- Read GWAS results ------")
          gwasResults <- private$read_gwas(pheno)

          message("------- Read LD  ------")
          ld <- private$read_ld(pheno)

          message("------- Read lead variants  ------")
          leadVar <- private$read_lead_variants(pheno)

          if (nrow(leadVar) == 0L) {
            message(sprintf("No lead variants found for %s. Skipping phenotype.", pheno))
            next
          }

          message("------- Filter GWAS results ------")
          # (kept as a message for parity; actual filtering happens per-locus below)

          # plot per lead variant (keep same loop + messages as original)
          for (i in seq_len(nrow(leadVar))) {
            # Per-locus guard so one bad locus doesn't stop the phenotype
            tryCatch({
              plot.name <- paste("RegionalPlot_", pheno, "_", leadVar$Locus_name[i], sep = "")
              # expected outputs (we’ll skip if any already exists)
              expected_png <- file.path(self$plot_dir, paste0(plot.name, ".png"))
              expected_pdf <- file.path(self$plot_dir, paste0(plot.name, ".pdf"))

              if (file.exists(expected_png) || file.exists(expected_pdf)) {
                message(sprintf("✓ Exists, skipping: %s (.png/.pdf)", plot.name))
                next
              }

              chr <- leadVar$CHROM[i]
              if (is.na(chr) || is.na(leadVar$GENPOS[i])) {
                message(sprintf("Lead variant row %d has NA CHROM/GENPOS; skipping %s.", i, plot.name))
                next
              }

              # subset GWAS to window
              plot.dat <- subset(
                gwasResults,
                CHROM == chr &
                  GENPOS > (leadVar$GENPOS[i] - self$plot_win) &
                  GENPOS < (leadVar$GENPOS[i] + self$plot_win)
              )

              if (nrow(plot.dat) == 0L) {
                message(sprintf("No GWAS points in window for %s; skipping.", plot.name))
                next
              }

              # add LD (merge on SNP; ld file is GRCh38-derived R2 but merged via "SNP" as in original)
              # be lenient if columns are missing
              if (!all(c("lead_snp", "SNP", "R2") %in% colnames(ld))) {
                message("LD file missing required columns (lead_snp/SNP/R2); setting R2 = 0 for all points.")
                plot.dat$R2 <- 0
              } else {
                ld_subset <- subset(ld, lead_snp == leadVar$SNP[i])[, c("SNP", "R2")]
                plot.dat <- merge(plot.dat, ld_subset, by = "SNP", all.x = TRUE)
                plot.dat$R2[is.na(plot.dat$R2)] <- 0
              }

              message(sprintf("Plotting: %s", plot.name))

              # rename columns to what makeRegionalPlot expects (as per original)
              names(plot.dat)[names(plot.dat) == "GENPOS"] <- "BP"
              names(plot.dat)[names(plot.dat) == "CHROM"]  <- "CHR"

              # Basic NA filtering to avoid plotting errors
              required_cols <- c("CHR", "BP", "P", "SNP")
              missing_req   <- setdiff(required_cols, colnames(plot.dat))
              if (length(missing_req) > 0L) {
                message(sprintf("Plot data missing required columns (%s); skipping %s.",
                                paste(missing_req, collapse = ", "), plot.name))
                next
              }
              plot.dat <- subset(plot.dat, !is.na(CHR) & !is.na(BP) & !is.na(P) & !is.na(SNP))
              if (nrow(plot.dat) == 0L) {
                message(sprintf("All points NA after filtering for %s; skipping.", plot.name))
                next
              }
              plot.dat <- as.data.frame(plot.dat)

              # recombination map: same behaviour — read file, standardise colnames
              recomb.dat <- tryCatch({
                rd <- read.table(self$genetic_map_path, sep = " ", h = TRUE)
                if (ncol(rd) < 4) stop("Genetic map must have 4 columns.")
                colnames(rd) <- c('Chromosome', 'Position(bp)', 'Rate(cM/Mb)', 'Map(cM)')
                as.data.frame(rd)
              }, error = function(e) {
                message("Warning: failed to read/parse genetic map: ", conditionMessage(e))
                NULL
              })

              # Call helper exactly as original if possible
              tryCatch({
                makeRegionalPlot(
                  plot.dat,
                  recomb.dat,
                  plot.name,
                  self$ensembl,
                  leadVar$SNP[i],
                  "SNP",
                  'P'
                )
              }, error = function(e) {
                message(sprintf("Error while plotting %s: %s. Skipping this plot.",
                                plot.name, conditionMessage(e)))
              })

            }, error = function(e) {
              message(sprintf("Unexpected error at locus index %d for phenotype %s: %s. Skipping this locus.",
                              i, pheno, conditionMessage(e)))
            })
          }

        }, error = function(e) {
          message(sprintf("Error while processing phenotype %s: %s. Skipping phenotype.",
                          pheno, conditionMessage(e)))
        })
      }

      message("Finished")
    }
  ),

  private = list(
    source_makeRegionalPlot = function() {
      # Source makeRegionalPlot.R from the working directory for portability
      makeRegionalPlot_path <- "makeRegionalPlot.R"
      if (file.exists(makeRegionalPlot_path)) {
        source(makeRegionalPlot_path)
      } else {
        stop("makeRegionalPlot.R not found in the working directory. Please ensure it is present for portability.")
      }
    },

    read_gwas = function(pheno) {
      path <- file.path(self$gwas_dir, paste0(pheno, "_38_37.txt"))
      tryCatch({
        dt <- fread(path)
        colnames(dt) <- make.names(colnames(dt))
        dt
      }, error = function(e) {
        stop(sprintf("Failed to read GWAS results for %s at %s: %s", pheno, path, conditionMessage(e)))
      })
    },

    read_ld = function(pheno) {
      tryCatch({
        dt <- as.data.frame(fread(self$ld_path, h = TRUE))
        dt <- dplyr::filter(dt, Phenotype == pheno)
        colnames(dt) <- make.names(colnames(dt))
        dt
      }, error = function(e) {
        message(sprintf("Warning: failed to read LD file at %s: %s. Proceeding with empty LD.",
                        self$ld_path, conditionMessage(e)))
        data.frame()
      })
    },

    read_lead_variants = function(pheno) {
      tryCatch({
        dt <- read.table(self$leadVar_path, h = TRUE, sep = ",")
        dplyr::filter(dt, Phenotype == pheno)
      }, error = function(e) {
        stop(sprintf("Failed to read lead variants at %s: %s", self$leadVar_path, conditionMessage(e)))
      })
    }
  )
)

# ----------------------------- Execute -----------------------------
gwas_dir        <- gwas_path
plot_dir        <- locuszoom_path
leadVar_dir     <- all_loci
ld_dir          <- var_file
genetic_map_dir <- genetic_map_hg38
phenotypes      <- traits

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

plotter <- LocusZoomPlotter$new(
  gwas_dir         = gwas_dir,
  plot_dir         = plot_dir,
  leadVar_path     = leadVar_dir,
  ld_path          = ld_dir,
  genetic_map_path = genetic_map_dir,
  phenotypes       = phenotypes,
  plot_win         = 500e3
)
plotter$run()
