# custom_coloc.R (OOP/R6 version; memory-safe & single-threaded)

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(coloc)
  library(data.table)
  library(stringr)
  library(R6)
})

# -------- Hard clamp threads & silence noisy systemd check --------
# (these take effect for most BLAS/LAPACK builds and data.table)
try({
  Sys.setenv(OMP_NUM_THREADS="1",
             OPENBLAS_NUM_THREADS="1",
             MKL_NUM_THREADS="1",
             NUMEXPR_NUM_THREADS="1",
             VECLIB_MAXIMUM_THREADS="1",
             R_DEFAULT_INTERNET_TIMEOUT="120")
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
  options(mc.cores = 1)
  # harmless on non-systemd images, prevents spam
  try(system("timedatectl 2>/dev/null", intern = TRUE), silent = TRUE)
}, silent = TRUE)

# Load pipeline config (expects: custom_gwas, all_loci, coloc_path, gwas_path, custom_trait_name)
source("config_R.R")

ColocRunner <- R6::R6Class(
  "ColocRunner",
  public = list(
    # config
    custom_gwas = NULL, all_loci = NULL, coloc_path = NULL, gwas_path = NULL, custom_trait_name = NULL,

    # data
    pheno_loci_dt = NULL,     # data.table for loci
    phenotypes    = NULL,
    custom_final  = NULL,     # data.table, preprocessed custom GWAS (min cols)

    initialize = function(custom_gwas, all_loci, coloc_path, gwas_path, custom_trait_name) {
      self$custom_gwas        <- custom_gwas
      self$all_loci           <- all_loci
      self$coloc_path         <- coloc_path
      self$gwas_path          <- gwas_path
      self$custom_trait_name  <- custom_trait_name

      if (!dir.exists(self$coloc_path)) {
        dir.create(self$coloc_path, recursive = TRUE, showWarnings = FALSE)
      }

      # read loci/lead SNPs for all phenotypes (CSV with header)
      # keep only the columns we actually use to reduce memory
      self$pheno_loci_dt <- fread(self$all_loci, showProgress = FALSE)

      # ensure required columns exist
      req_cols <- c("Phenotype","SNP","CHROM","GENPOS")
      missing <- setdiff(req_cols, names(self$pheno_loci_dt))
      if (length(missing)) stop("Missing required columns in all_loci: ", paste(missing, collapse=", "))

      # normalize types
      self$pheno_loci_dt[, CHROM := as.character(CHROM)]
      self$pheno_loci_dt[, GENPOS := as.integer(GENPOS)]

      self$phenotypes <- unique(self$pheno_loci_dt$Phenotype)

      # read & prepare custom GWAS once (select minimal columns)
      self$custom_final  <- self$prepare_custom_gwas_(self$custom_gwas)
    },

    run = function() {
      for (p in self$phenotypes) {
        cat(str_interp('Processing ${p}\n'))
        tryCatch({
          self$run_one_phenotype_(p)
          gc(verbose = FALSE)
        }, error = function(e) {
          err_dir <- file.path(self$coloc_path, "Errors")
          if (!dir.exists(err_dir)) dir.create(err_dir, recursive = TRUE, showWarnings = FALSE)
          cat(str_interp("Error while running coloc for ${p}: "), e$message, "\n")
          saveRDS(e$message, file.path(err_dir, str_interp('error_message_${self$custom_trait_name}_${p}.rds')))
        })
      }
      invisible(TRUE)
    },

    # -------------------------- internals --------------------------

    prepare_custom_gwas_ = function(path) {
      # Minimal columns required downstream
      # chromosome, base_pair_location, effect_allele, other_allele,
      # p_value, beta, standard_error, effect_allele_frequency, base_pair_location37 (optional)
      sel <- c("chromosome","base_pair_location","effect_allele","other_allele",
               "p_value","beta","standard_error","effect_allele_frequency","base_pair_location37")
      dt <- suppressWarnings(fread(path, select = intersect(sel, names(fread(path, nrows = 0))),
                                   showProgress = FALSE))

      # type coercions
      setnames(dt, old = intersect(c("beta"), names(dt)), new = intersect(c("beta"), names(dt)))
      if ("base_pair_location" %in% names(dt)) dt[, base_pair_location := as.integer(base_pair_location)]
      if ("base_pair_location37" %in% names(dt)) dt[, base_pair_location37 := as.integer(base_pair_location37)]
      if ("effect_allele_frequency" %in% names(dt)) {
        dt[, effect_allele_frequency := as.numeric(effect_allele_frequency)]
      }

      # build IDs (flip-aware) with data.table
      dt[, newID1 := paste(chromosome, base_pair_location, effect_allele, other_allele, sep=":")]
      dt[, newID2 := paste(chromosome, base_pair_location, other_allele, effect_allele, sep=":")]

      # keep one row per newID1 (smallest p_value) – do this with data.table for memory efficiency
      if ("p_value" %in% names(dt)) {
        setorder(dt, p_value)
      }
      dt <- dt[!duplicated(newID1)]

      # ensure no missing positions
      dt <- dt[!is.na(base_pair_location)]

      # keys for fast join later
      setkeyv(dt, "newID1")
      dt
    },

    run_one_phenotype_ = function(p) {
      # read phenotype GWAS for this phenotype only; prune columns to those used
      gwas_file <- file.path(self$gwas_path, str_interp('${p}_38_37.txt'))
      sel <- c("CHROM","GENPOS","BETA","SE","SNP","ALLELE0","ALLELE1","A1FREQ","INFO","MAF","N")
      gwas_data <- fread(gwas_file, select = intersect(sel, names(fread(gwas_file, nrows = 0))),
                         showProgress = FALSE)
      gwas_data[, CHROM := as.character(CHROM)]
      gwas_data[, GENPOS := as.integer(GENPOS)]

      # 1MB blocks around each lead SNP
      lead_vars <- self$pheno_loci_dt[Phenotype == p, unique(SNP)]
      gwas_data_loci <- self$pheno_loci_dt[Phenotype == p, .(SNP, CHROM, GENPOS)]

      gwas_blocks <- self$build_blocks_(lead_vars, gwas_data_loci, gwas_data)
      rm(gwas_data); gc(FALSE)

      gwas_common <- self$join_blocks_with_custom_(lead_vars, gwas_blocks, self$custom_final)
      rm(gwas_blocks); gc(FALSE)

      # coloc per lead variant
      coloc_res <- self$run_coloc_(lead_vars, gwas_common)

      # filter H4/H3
      coloc_h4_ge08 <- Filter(function(x) x$summary[[6]] > 0.8, coloc_res)
      coloc_h3_ge08 <- Filter(function(x) x$summary[[5]] > 0.8, coloc_res)

      out_path_yes  <- file.path(self$coloc_path, str_interp('coloc_${self$custom_trait_name}_res_all_${p}.rds'))
      out_path_no   <- file.path(self$coloc_path, str_interp('coloc_${self$custom_trait_name}_res_all_noPP4${p}.rds'))

      if (length(coloc_h4_ge08) >= 1) {
        coloc_res_dat <- self$assemble_pp4_dat_(p, coloc_h4_ge08, gwas_common)
        saveRDS(list(Full_res = coloc_res, H4 = coloc_h4_ge08, H3 = coloc_h3_ge08, PP4 = coloc_res_dat), out_path_yes)
      } else {
        saveRDS(list(Full_res = coloc_res, H4 = coloc_h4_ge08, H3 = coloc_h3_ge08), out_path_no)
      }

      invisible(TRUE)
    },

    build_blocks_ = function(lead_vars, gwas_data_loci, gwas_data) {
      # Prepare newID1 in gwas_data once (saves repeated paste)
      gwas_data[, newID1 := paste(CHROM, GENPOS, ALLELE0, ALLELE1, sep=":")]
      setkeyv(gwas_data, c("CHROM","GENPOS"))

      # For each lead SNP get +/- 500kb block on matching chromosome
      res_list <- vector("list", length(lead_vars)); names(res_list) <- lead_vars
      for (k in seq_along(lead_vars)) {
        snp <- lead_vars[[k]]
        chr <- gwas_data_loci[SNP == snp, CHROM][1]
        pos <- gwas_data_loci[SNP == snp, GENPOS][1]
        blk <- gwas_data[CHROM == chr & GENPOS >= (pos - 5e5) & GENPOS <= (pos + 5e5)]
        res_list[[k]] <- blk
      }
      res_list
    },

    join_blocks_with_custom_ = function(lead_vars, gwas_blocks, custom_final) {
      out <- vector("list", length(lead_vars)); names(out) <- lead_vars

      for (k in seq_along(lead_vars)) {
        i <- lead_vars[[k]]
        cat(str_interp('Processing lead variant ${i}\n'))
        block <- copy(gwas_blocks[[k]])

        # DT joins: by newID1 and flipped newID2
        setkeyv(block, "newID1")
        res1 <- custom_final[block, on=.(newID1), nomatch = 0L]      # join where custom$newID1 == block$newID1
        # for flipped: match custom$newID2 to block$newID1
        setkeyv(custom_final, "newID2")
        res2 <- custom_final[block, on=.(newID2 = newID1), nomatch = 0L]
        setkeyv(custom_final, "newID1")  # restore

        # phenotype (quant)
        # Keep only once per snp later in run_coloc_
        phenotype_gwas <- rbindlist(list(res1, res2), use.names = TRUE, fill = TRUE)
        phenotype_gwas[, `:=`(
          varbeta = SE^2,
          position = as.integer(GENPOS),
          type = 'quant',
          PHENOTYPE_MAF = MAF,
          PHENOTYPE_p = P,
          PHENOTYPE_SE = SE,
          PHENOTYPE_N = N,
          Org_ID = SNP,
          PHENOTYPE_NEA = ALLELE0,
          PHENOTYPE_EA  = ALLELE1,
          PHENOTYPE_EAF = A1FREQ,
          PHENOTYPE_INFO = INFO
        )]
        setnames(phenotype_gwas, old = "newID1", new = "snp")
        phenotype_gwas <- phenotype_gwas[, .(beta = BETA, varbeta, snp, position, type,
                                             PHENOTYPE_MAF, PHENOTYPE_p, PHENOTYPE_SE,
                                             PHENOTYPE_N, Org_ID, PHENOTYPE_NEA,
                                             PHENOTYPE_EA, PHENOTYPE_EAF, PHENOTYPE_INFO)]

        # custom (cc), flip for matches from res1
        # res1 are cases where custom alleles == block alleles; to align with phenotype,
        # we flip custom beta and allele frequency (as in your original)
        if (nrow(res1)) {
          res1[, `:=`(beta = -beta, effect_allele_frequency = 1 - effect_allele_frequency)]
        }
        cust <- rbindlist(list(res1, res2), use.names = TRUE, fill = TRUE)
        cust[, `:=`(
          MAF = pmin(effect_allele_frequency, 1 - effect_allele_frequency),
          varbeta = (standard_error)^2,
          position = as.integer(base_pair_location),
          type = 'cc'
        )]
        setnames(cust, old = "newID1", new = "snp")
        custom <- cust[, .(beta, varbeta, snp, position, type,
                           CUSTOM_MAF = MAF, CUSTOM_SE = standard_error, CUSTOM_P = p_value)]
        # optional GRCh37 passthrough (kept as before): attach if present in source
        if ("base_pair_location37" %in% names(cust)) {
          custom[, position37 := as.integer(base_pair_location37)]
        }

        out[[k]] <- list(phenotype_gwas, custom)
        rm(res1, res2, cust, custom, phenotype_gwas, block); gc(FALSE)
      }
      out
    },

    run_coloc_ = function(lead_vars, gwas_common) {
      # single-threaded; deduplicate snps before coloc
      res <- vector("list", length(lead_vars)); names(res) <- lead_vars
      for (k in seq_along(lead_vars)) {
        i <- lead_vars[[k]]

        D1_df <- as.data.table(gwas_common[[i]][[1]])
        setkeyv(D1_df, "snp"); D1_df <- unique(D1_df, by="snp"); setorder(D1_df, snp)
        D1 <- list(D1_df$beta, D1_df$varbeta, D1_df$snp, D1_df$position, 'quant', 1)
        names(D1) <- c('beta','varbeta','snp','position','type','sdY')

        D2_df <- as.data.table(gwas_common[[i]][[2]])
        setkeyv(D2_df, "snp"); D2_df <- unique(D2_df, by="snp"); setorder(D2_df, snp)
        D2 <- list(D2_df$beta, D2_df$varbeta, D2_df$snp, D2_df$position, 'cc')
        names(D2) <- c('beta','varbeta','snp','position','type')

        res[[k]] <- coloc.abf(dataset1 = D1, dataset2 = D2)
        rm(D1_df, D2_df, D1, D2); gc(FALSE)
      }
      res
    },

    assemble_pp4_dat_ = function(phenotype_id, coloc_h4_ge08, gwas_common) {
      out <- vector("list", length(coloc_h4_ge08)); names(out) <- names(coloc_h4_ge08)
      for (j in seq_along(out)) {
        n <- names(out)[j]
        left1 <- self$pheno_loci_dt[Phenotype == phenotype_id & SNP == n,
                                    .(Phenotype, Locus_name, Lead_variant = SNP)]
        res_tbl <- coloc_h4_ge08[[n]]$results
        # safe joins with data.table
        D1 <- as.data.table(gwas_common[[n]][[1]])
        D2 <- as.data.table(gwas_common[[n]][[2]])
        setkeyv(D1, "snp"); setkeyv(D2, "snp")
        res_dt <- as.data.table(res_tbl)
        setkeyv(res_dt, "snp")
        merged <- res_dt[D1, nomatch=0][D2, nomatch=0]
        out[[j]] <- cbind(left1, merged)
        rm(D1, D2, res_dt, merged); gc(FALSE)
      }
      out
    }
  )
)

# ------------------------------- Execute -------------------------------
runner <- ColocRunner$new(
  custom_gwas        = custom_gwas,
  all_loci           = all_loci,
  coloc_path         = coloc_path,
  gwas_path          = gwas_path,
  custom_trait_name  = custom_trait_name
)
invisible(runner$run())
