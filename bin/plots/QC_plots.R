#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(R6)
  library(tidyverse)
  library(magrittr)
  library(data.table)
  library(patchwork)
  library(here)
})

source("config_R.R")

GWASPlotter <- R6::R6Class(
  "GWASPlotter",
  public = list(
    gwas_path = NULL,
    plot_outpath = NULL,
    cutoff = 5e-08,
    qq_dir = NULL,
    qc_dir = NULL,

    initialize = function(gwas_path, plot_outpath, cutoff = 5e-08) {
      self$gwas_path    <- gwas_path
      self$plot_outpath <- plot_outpath
      self$cutoff       <- cutoff

      self$qq_dir <- file.path(self$plot_outpath, "QQ_plots")
      self$qc_dir <- file.path(self$plot_outpath, "QC_plots")

      private$make_dir(self$qq_dir)
      private$make_dir(self$qc_dir)
    },

    read_trait = function(trait) {
      file_name <- paste0(trait, "_38_37_rsids.txt")
      full_path <- file.path(self$gwas_path, file_name)
      message("[", trait, "] Reading: ", full_path)
      data.table::fread(
        full_path,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        fill = TRUE
      )
    },

    qqplot_save = function(res_data, pheno_name) {
      gg_qqplot <- function(ps, ci = 0.95) {
        n  <- length(ps)
        df <- data.frame(
          observed = -log10(sort(ps)),
          expected = -log10(ppoints(n)),
          clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
          cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
        )
        log10Pe <- expression(paste("Expected -log"[10], plain(P)))
        log10Po <- expression(paste("Observed -log"[10], plain(P)))
        ggplot(df) +
          geom_ribbon(aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1) +
          geom_point(aes(expected, observed), shape = 1, size = 3) +
          geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
          xlab(log10Pe) + ylab(log10Po)
      }

      inflation <- function(ps) {
        chisq <- qchisq(1 - ps, 1)
        median(chisq) / qchisq(0.5, 1)
      }

      qqplot_obj <- gg_qqplot(res_data[["P"]]) +
        theme_bw(base_size = 20) +
        annotate(
          geom = "text",
          x = -Inf, y = Inf, hjust = -0.15, vjust = 1 + 0.45,
          label = sprintf("λ = %.3f", inflation(res_data[["P"]])),
          size = 6
        ) +
        theme(axis.ticks = element_line(size = 0.5), panel.grid = element_blank()) +
        ggtitle(pheno_name)

      # Save ONLY into QQ_plots using the `path=` argument
      ggplot2::ggsave(
        plot = qqplot_obj,
        filename = stringr::str_interp('QQ_plot_${pheno_name}.png'),
        path = self$qq_dir,
        width = 10, height = 8, dpi = 300
      )
      rm(qqplot_obj)
      invisible(NULL)
    },

    qc_plots_save = function(res_data, pheno_name) {
      plot1 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(MAF, `-log10 P`, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() +
        theme(legend.position = 'none',
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12)) +
        ggtitle(pheno_name)

      plot2 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(MAF, BETA, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() + theme(legend.position = 'none')

      plot3 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(MAF, SE, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() + theme(legend.position = 'none')

      plot4 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(BETA, `-log10 P`, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() + theme(legend.position = 'none')

      plot5 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(BETA, SE, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() + theme(legend.position = 'none')

      plot6 <- ggplot(res_data %>% mutate(`-log10 P` = -log10(P)) %>%
                        mutate(MAF_cutoff = ifelse(MAF >= 0.05, 'Yes', 'No'),
                               Significant = ifelse(P < self$cutoff, 'Yes', 'No')),
                      aes(SE, `-log10 P`, colour = MAF_cutoff, shape = Significant)) +
        geom_point(size = 2, alpha = 0.8) +
        theme_classic() + theme(legend.position = 'none')

      plot_final <- plot1 + plot2 + plot3 + plot4 + plot5 + plot6 + plot_layout(ncol = 2)

      # Save ONLY into QC_plots using the `path=` argument
      ggplot2::ggsave(
        plot = plot_final,
        filename = paste0(pheno_name, '_QC_plots.png'),
        path = self$qc_dir,
        width = 15, height = 10, dpi = 300
      )
      rm(plot_final)
      invisible(NULL)
    },

    process_trait = function(trait) {
      res <- self$read_trait(trait)
      print(trait)
      self$qqplot_save(res, trait)
      print('qq plot finished')
      self$qc_plots_save(res, trait)
      print('qc plot finished')
      invisible(NULL)
    },

    run = function(traits) {
      purrr::walk(traits, self$process_trait)
      invisible(NULL)
    }
  ),

  private = list(
    make_dir = function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  )
)

# ----------------------------- Execute -----------------------------
plotter <- GWASPlotter$new(gwas_path = gwas_path, plot_outpath = plot_outpath, cutoff = 5e-08)
plotter$run(traits)
