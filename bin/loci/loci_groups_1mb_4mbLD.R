#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(igraph)
  library(here)
  library(R6)
})

# We group SNPs into a locus if they’re within ±500 kb or if they’re in LD R²≥0.4 within 4 Mb, on the same chromosome. We build one union graph of these edges and take connected components as loci. This guarantees that ±500 kb neighbours always group, and LD can merge separate distance clusters (bounded to 4 Mb), with transitivity handled correctly.


source("config_R.R")

LocusGrouper <- R6::R6Class(
  "LocusGrouper",
  public = list(
    all_loci    = NULL,
    var_file    = NULL,
    traits      = NULL,
    output_path = NULL,

    # params
    DIST_MAX_BP = 5e5,     # distance rule: <= 500 kb
    LD_R2_MIN   = 0.4,     # LD rule threshold
    LD_MAX_BP   = 4e6,     # allow LD edges up to 4 Mb

    # data
    df_all   = NULL,
    snp_pos  = NULL,   # Lead_SNP, CHROM, GENPOS, KeyPos
    ld       = NULL,
    snpA_col = NULL,
    snpB_col = NULL,

    edges_dist  = NULL,
    edges_ld    = NULL,
    edges_union = NULL,

    g          = NULL,
    locus_map  = NULL,
    output     = NULL,

    initialize = function(all_loci, var_file, traits, output_path) {
      self$all_loci    <- all_loci
      self$var_file    <- var_file
      self$traits      <- traits
      self$output_path <- output_path
    },

    run = function() {
      self$load_inputs()
      self$prepare_positions()
      self$load_ld()
      self$detect_snp_cols()
      self$build_distance_edges()
      self$build_ld_edges_with_window()
      self$union_edges_and_components()
      self$attach_reduce_and_format()
      self$normalize_order_and_renumber()
      self$write_output()
    },

    load_inputs = function() {
      self$df_all <- fread(self$all_loci)
      stopifnot("SNP" %in% names(self$df_all))
      setnames(self$df_all, "SNP", "Lead_SNP")
    },

    prepare_positions = function() {
      pos_cols <- c("Lead_SNP", "CHROM", "GENPOS")
      self$snp_pos <- unique(self$df_all[, ..pos_cols])
      self$snp_pos[, CHROM := as.numeric(CHROM)]
      self$snp_pos[, GENPOS := as.numeric(GENPOS)]
      self$snp_pos[, KeyPos := paste0(CHROM, ":", GENPOS)]
      # keys for foverlaps
      self$snp_pos[, `:=`(start = GENPOS, end = GENPOS)]
      setkey(self$snp_pos, CHROM, start, end)
    },

    load_ld = function() {
      self$ld <- fread(self$var_file)
      if ("Gene_Symbol" %in% names(self$ld)) {
        self$ld <- self$ld[Gene_Symbol != ""]
      }
      if (!("R2" %in% names(self$ld))) stop("LD file must contain an 'R2' column.")
    },

    detect_snp_cols = function() {
      cand_pairs <- list(
        c("SNP_A","SNP_B"),
        c("SNP1","SNP2"),
        c("SNP_A1","SNP_A2"),
        c("VAR1","VAR2"),
        c("ID1","ID2"),
        c("ID_A","ID_B")
      )
      found <- NULL
      for (p in cand_pairs) if (all(p %in% names(self$ld))) { found <- p; break }
      if (is.null(found)) {
        if (ncol(self$ld) < 3) stop("LD file must have at least 3 columns (two SNP ids + R2).")
        found <- names(self$ld)[1:2]
      }
      self$snpA_col <- found[1]
      self$snpB_col <- found[2]
    },

    # Build edges from the distance rule (<= 500 kb on same chromosome)
    build_distance_edges = function() {
      # Build ± window intervals around each position (table A)
      a <- self$snp_pos[, .(Lead_SNP, CHROM,
                             start = pmax(GENPOS - self$DIST_MAX_BP, 0),
                             end   = GENPOS + self$DIST_MAX_BP)]
      setkey(a, CHROM, start, end)

      # Points table at exact positions (table B)
      b <- self$snp_pos[, .(Lead_SNP, CHROM, start, end)]
      setkey(b, CHROM, start, end)

      # NOTE: some data.table versions don't support allow.cartesian in foverlaps
      ov <- foverlaps(b, a, nomatch = 0L)

      # unique undirected edges; drop self-loops
      edges <- ov[Lead_SNP != i.Lead_SNP,
                  .(from = pmin(Lead_SNP, i.Lead_SNP),
                    to   = pmax(Lead_SNP, i.Lead_SNP))]
      self$edges_dist <- unique(edges)
    },

    .norm_key = function(x) {
      x <- as.character(x)
      x <- gsub("^chr", "", x, ignore.case = TRUE)
      sub("^([^:]+:[0-9]+).*$", "\\1", x)
    },

    # Build edges from LD rule: same chr, r2 >= 0.4, and distance <= 4 Mb
    build_ld_edges_with_window = function() {
      snpA <- self$snpA_col; snpB <- self$snpB_col
      ld_f <- self$ld[R2 >= self$LD_R2_MIN, .SD, .SDcols = c(snpA, snpB, "R2")]
      if (nrow(ld_f) == 0L) { self$edges_ld <- data.table(from=character(), to=character()); return(invisible(NULL)) }

      setDT(ld_f)
      ld_f[, KeyA := self$.norm_key(get(snpA))]
      ld_f[, KeyB := self$.norm_key(get(snpB))]

      # Only keep pairs we can position-match
      ld_f <- ld_f[KeyA %in% self$snp_pos$KeyPos & KeyB %in% self$snp_pos$KeyPos]
      if (nrow(ld_f) == 0L) { self$edges_ld <- data.table(from=character(), to=character()); return(invisible(NULL)) }

      mapA <- self$snp_pos[, .(KeyA = KeyPos, Lead_SNP_A = Lead_SNP, CHROM_A = CHROM, GENPOS_A = GENPOS)]
      mapB <- self$snp_pos[, .(KeyB = KeyPos, Lead_SNP_B = Lead_SNP, CHROM_B = CHROM, GENPOS_B = GENPOS)]

      ld_f <- merge(ld_f, mapA, by = "KeyA", allow.cartesian = TRUE)
      ld_f <- merge(ld_f, mapB, by = "KeyB", allow.cartesian = TRUE)

      # same chromosome and within LD_MAX_BP
      ld_f <- ld_f[CHROM_A == CHROM_B & abs(GENPOS_A - GENPOS_B) <= self$LD_MAX_BP]

      if (nrow(ld_f) == 0L) { self$edges_ld <- data.table(from=character(), to=character()); return(invisible(NULL)) }

      edges <- ld_f[, .(from = pmin(Lead_SNP_A, Lead_SNP_B),
                        to   = pmax(Lead_SNP_A, Lead_SNP_B))]
      edges <- unique(edges[from != to])
      self$edges_ld <- edges
    },

    union_edges_and_components = function() {
      all_vertices <- data.frame(name = unique(self$df_all$Lead_SNP), stringsAsFactors = FALSE)

      # union of distance AND LD edges (deduplicated)
      edges <- rbindlist(list(self$edges_dist, self$edges_ld), use.names = TRUE, fill = TRUE)
      edges <- unique(edges[complete.cases(edges)])

      if (nrow(edges) > 0L) {
        self$g <- graph_from_data_frame(d = edges, vertices = all_vertices, directed = FALSE)
      } else {
        self$g <- make_empty_graph(n = nrow(all_vertices), directed = FALSE)
        self$g <- set_vertex_attr(self$g, "name", value = all_vertices$name)
      }

      memb <- components(self$g)$membership
      self$locus_map <- data.table(Lead_SNP = names(memb), Locus_number = as.integer(memb))
    },

    attach_reduce_and_format = function() {
      out <- merge(self$df_all, self$locus_map, by = "Lead_SNP", all.x = TRUE)
      out <- out %>% dplyr::relocate(Locus_number, .after = Lead_SNP)
      setDT(out)
      # one row per (Phenotype, Lead_SNP), keeping the smallest P
      out <- out[order(P), .SD[1], by = .(Phenotype, Lead_SNP)]
      out[, Phenotype_Count_per_Locus_name := uniqueN(Phenotype), by = Locus_name]

      out <- dplyr::select(
        out,
        Phenotype, Phenotype_Count_per_Locus_name, Locus_name, Locus_number,
        Lead_SNP, CHROM, GENPOS, GENPOS_hg19, ALLELE0, ALLELE1, MAF, BETA, SE, P, Method
      )
      self$output <- out
    },

    normalize_order_and_renumber = function() {
      phenotype_order <- stringr::str_squish(toupper(self$traits))
      self$output[, Phenotype := as.character(Phenotype)]
      self$output[, Phenotype := stringr::str_squish(Phenotype)]
      self$output[, Phenotype := toupper(Phenotype)]
      is_match <- self$output$Phenotype %in% phenotype_order
      self$output[is_match, Phenotype := factor(Phenotype, levels = phenotype_order)]
      self$output <- self$output[order(Phenotype, as.numeric(CHROM), GENPOS)]

      # compact 1..K
      self$output[, Locus_number := .GRP, by = Locus_number]
    },

    write_output = function() {
      fwrite(self$output, self$output_path)
    }
  )
)

# ----------------------------- Execute -----------------------------
grouper <- LocusGrouper$new(
  all_loci = all_loci,
  var_file = var_file,
  traits = traits,
  output_path = locus_blocks_output
)

grouper$run()
