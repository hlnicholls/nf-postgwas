# process_custom_coloc_output.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(vautils)
  library(stringr)
  library(data.table)
  library(rlang)
  library(R6)
})

source("config_R.R")

# ----------------------------- R6: ColocProcessor -----------------------------

ColocProcessor <- R6::R6Class(
  "ColocProcessor",
  public = list(
    coloc_path = NULL,
    all_loci_path = NULL,
    var_file = NULL,
    custom_trait_name = NULL,
    file_list = NULL,
    phenotypes = NULL,
    pheno_loci_df = NULL,
    custom_pattern = NULL,
    ct_beta = NULL,
    ct_se = NULL,
    ct_p = NULL,

    initialize = function(coloc_path, all_loci_path, var_file, custom_trait_name) {
      self$coloc_path        <- normalizePath(coloc_path, mustWork = FALSE)
      self$all_loci_path     <- all_loci_path
      self$var_file          <- var_file
      self$custom_trait_name <- custom_trait_name

      if (!dir.exists(self$coloc_path)) {
        dir.create(self$coloc_path, recursive = TRUE, showWarnings = FALSE)
      }

      self$file_list <- list.files(self$coloc_path, pattern = "\\.rds$", full.names = TRUE)
      self$pheno_loci_df <- data.table::fread(self$all_loci_path)
      self$phenotypes <- unique(self$pheno_loci_df$Phenotype)

      self$custom_pattern <- paste0("coloc_", self$custom_trait_name, "_res_all_")
      self$ct_beta <- paste0(self$custom_trait_name, "_BETA")
      self$ct_se   <- paste0(self$custom_trait_name, "_SE")
      self$ct_p    <- paste0(self$custom_trait_name, "_P")
    },

    run = function() {
      colocalised <- self$get_colocalised_phenotypes()

      if (length(colocalised) == 0) {
        cat('No colocalised phenotypes found. Creating empty output file and exiting.\n')
        data.table::fwrite(data.frame(), file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_all_variants_pp4.csv')))
        data.table::fwrite(data.frame(), file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_pp4_all.csv')))
        data.table::fwrite(data.frame(), file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_pp3_all.csv')))
        return(invisible(TRUE))
      }

      # Build per-locus top-SNP PP4 table and write
      pp4_top <- self$build_pp4_top(colocalised) %>%
        self$attach_ld_snps()
      data.table::fwrite(pp4_top, file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_pp4_all.csv')))

      # Build all-variants PP4 table and write
      pp4_all <- self$build_pp4_all_variants()
      data.table::fwrite(pp4_all, file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_all_variants_pp4.csv')))

      # Build PP3 (H3) table and write
      pp3_all <- self$build_pp3_all()
      data.table::fwrite(pp3_all, file.path(self$coloc_path, paste0(self$custom_trait_name, '_coloc_pp3_all.csv')))

      invisible(TRUE)
    },

    # ----------------------------- Builders -----------------------------

    build_pp4_top = function(colocalised) {
      out <- map(colocalised, function(p) {
        trait <- readRDS(file.path(self$coloc_path, str_interp('coloc_${self$custom_trait_name}_res_all_${p}.rds')))

        # All colocalised rows (unused downstream, but keep for parity with original code flow)
        trait_PP4_all <- map(trait$PP4, function(x) x %>%
                               arrange(desc(SNP.PP.H4))) %>%
          bind_rows() %>%
          dplyr::select(
            Phenotype, Locus_name, Lead_variant,
            Variant = snp, GRCh38_pos = position.x,
            EA = PHENOTYPE_EA, NEA = PHENOTYPE_NEA,
            PHENOTYPE_EAF, PHENOTYPE_INFO,
            PHENOTYPE_BETA = beta.x, PHENOTYPE_SE, PHENOTYPE_P = PHENOTYPE_p,
            CUSTOM_BETA = beta.y, CUSTOM_SE, CUSTOM_P,
            V.df1, z.df1, r.df1, lABF.df1, V.df2, z.df2, r.df2, lABF.df2,
            internal.sum.lABF, SNP.PP.H4
          )

        # Top SNP per locus
        trait_PP4_top <- map(trait$PP4, function(x) x %>%
                               arrange(desc(SNP.PP.H4)) %>%
                               slice(1)) %>%
          bind_rows() %>%
          dplyr::select(
            Phenotype, Locus_name, Lead_variant,
            Variant = snp, GRCh38_pos = position.x,
            EA = PHENOTYPE_EA, NEA = PHENOTYPE_NEA,
            PHENOTYPE_EAF, PHENOTYPE_INFO,
            PHENOTYPE_BETA = beta.x, PHENOTYPE_SE, PHENOTYPE_P = PHENOTYPE_p,
            CUSTOM_BETA = beta.y, CUSTOM_SE, CUSTOM_P,
            V.df1, z.df1, r.df1, lABF.df1, V.df2, z.df2, r.df2, lABF.df2,
            internal.sum.lABF, SNP.PP.H4
          )

        # Optional GRCh37 augmentation
        trait_PP4_top <- self$maybe_attach_grch37(trait_PP4_top, trait)

        # Rename CUSTOM_* → <trait>_*
        trait_PP4_top <- trait_PP4_top %>%
          dplyr::rename(
            !!sym(self$ct_beta) := CUSTOM_BETA,
            !!sym(self$ct_se)   := CUSTOM_SE,
            !!sym(self$ct_p)    := CUSTOM_P
          )

        trait_PP4_top
      }) %>% bind_rows()

      out %>% dplyr::rename(SNP = Variant, lead_snp = Lead_variant)
    },

    build_pp4_all_variants = function() {
      rds_files <- list.files(self$coloc_path, pattern = "\\.rds$", full.names = TRUE)
      rds_files <- rds_files[grepl(self$custom_pattern, rds_files) & !grepl("noPP4", rds_files)]

      all_PP4_SNPs <- data.frame()

      for (file in rds_files) {
        trait <- readRDS(file)

        trait_PP4 <- map(trait$PP4, function(x) x %>%
                           arrange(desc(SNP.PP.H4))) %>%
          bind_rows() %>%
          dplyr::select(
            Phenotype, Locus_name, Lead_variant,
            Variant = snp, GRCh38_pos = position.x,
            EA = PHENOTYPE_EA, NEA = PHENOTYPE_NEA,
            PHENOTYPE_EAF, PHENOTYPE_INFO,
            PHENOTYPE_BETA = beta.x, PHENOTYPE_SE, PHENOTYPE_P = PHENOTYPE_p,
            CUSTOM_BETA = beta.y, CUSTOM_SE, CUSTOM_P,
            V.df1, z.df1, r.df1, lABF.df1, V.df2, z.df2, r.df2, lABF.df2,
            internal.sum.lABF, SNP.PP.H4
          )

        trait_PP4 <- self$maybe_attach_grch37(trait_PP4, trait)

        trait_PP4 <- trait_PP4 %>%
          dplyr::rename(
            !!sym(self$ct_beta) := CUSTOM_BETA,
            !!sym(self$ct_se)   := CUSTOM_SE,
            !!sym(self$ct_p)    := CUSTOM_P
          ) %>%
          mutate(CHROM = sapply(strsplit(as.character(Variant), ":"), function(x) x[1]))

        geneN <- find_nearest_gene(trait_PP4, flanking = 10, build = "hg38",
                                   collapse = FALSE, snp = "Variant", chr = "CHROM", bp = "GRCh38_pos") %>%
          dplyr::rename(Variant = rsid) %>%
          dplyr::mutate(distance = ifelse(distance == "intergenic", 0, distance))

        trait_PP4 <- merge(
          subset(geneN[order(geneN$Variant, abs(as.numeric(geneN$distance))), ],
                 !duplicated(Variant))[, c("Variant", "GENE")],
          trait_PP4, by = "Variant"
        ) %>%
          dplyr::select(Nearest_Gene_10kb = GENE, dplyr::everything())

        all_PP4_SNPs <- rbind(all_PP4_SNPs, trait_PP4)
      }

      all_PP4_SNPs
    },

    build_pp3_all = function() {
      rds_files <- list.files(self$coloc_path, pattern = "\\.rds$", full.names = TRUE)
      rds_files <- rds_files[grepl(self$custom_pattern, rds_files)]

      H3_SNPs <- data.frame(Variant = character(), PP.H3.abf = numeric(), Phenotype = character(), stringsAsFactors = FALSE)

      for (file in rds_files) {
        phenotype <- stringr::str_extract(basename(file), paste0("(?<=", self$custom_pattern, ").*(?=\\.rds)"))
        trait <- readRDS(file)

        if (is.list(trait$H3)) {
          for (snp in names(trait$H3)) {
            if (is.list(trait$H3[[snp]]) &&
                !is.null(trait$H3[[snp]]$summary) &&
                "PP.H3.abf" %in% names(trait$H3[[snp]]$summary)) {

              pp_h3_abf <- trait$H3[[snp]]$summary["PP.H3.abf"]
              H3_SNPs <- rbind(
                H3_SNPs,
                data.frame(Variant = snp, PP.H3.abf = pp_h3_abf, Phenotype = phenotype, stringsAsFactors = FALSE)
              )
            }
          }
        }
      }

      H3_SNPs
    },

    # ----------------------------- Helpers -----------------------------

    get_colocalised_phenotypes = function() {
      colocalised <- self$file_list[grepl(self$custom_pattern, self$file_list) & !grepl('noPP4', self$file_list)] %>%
        basename() %>%
        gsub(paste0(self$custom_pattern, '|.rds'), '', .) %>%
        unique()

      colocalised <- colocalised[order(match(colocalised, self$phenotypes))]
      colocalised[colocalised != ""]
    },

    maybe_attach_grch37 = function(df, trait_obj) {
      if (!length(trait_obj$PP4)) return(df)

      # Check the first element for the presence of position37
      first_has_37 <- isTRUE("position37" %in% names(trait_obj$PP4[[1]]))
      if (!first_has_37) return(df)

      df <- df %>% mutate(GRCh37_pos = NA_integer_)
      for (i in seq_along(trait_obj$PP4)) {
        locus_data <- trait_obj$PP4[[i]]
        if ("position37" %in% names(locus_data)) {
          # Use Locus_name to align the block (as in the original script)
          df$GRCh37_pos[df$Locus_name == locus_data$Locus_name[1]] <- locus_data$position37
        }
      }
      df %>% dplyr::select(Phenotype, Locus_name, Lead_variant, Variant, GRCh38_pos, GRCh37_pos, dplyr::everything())
    },

    attach_ld_snps = function(pp4_df) {
      ld_snps <- data.table::fread(self$var_file) %>%
        dplyr::select(SNP, Gene_Symbol, rsid_1kg)

      merge(pp4_df, ld_snps, by = "SNP", all.x = TRUE) %>%
        dplyr::rename(Nearest_Gene_10kb = Gene_Symbol)
    }
  )
)

# ----------------------------- Execute -----------------------------

processor <- ColocProcessor$new(
  coloc_path        = coloc_path,
  all_loci_path     = all_loci,
  var_file          = var_file,
  custom_trait_name = custom_trait_name
)

invisible(processor$run())
