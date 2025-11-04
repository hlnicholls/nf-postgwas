#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(here)
  library(stringr)
})

# -------- SETTINGS --------
genomewide <- 5e-8
suggestive <- 1e-5
label_dist_kb <- 500
save_width_px <- 4500
save_height_px <- 2500
save_res <- 300
autosomes_only <- TRUE
options(datatable.integer64 = "numeric")
# --------------------------

source("config_R.R")
dir.create(plot_outpath, showWarnings = FALSE, recursive = TRUE)

# Palette fixed across chr 1..22
chr_levels <- as.character(1:22)
chr_palette <- setNames(
  hcl(h = seq(15, 375, length.out = 23)[-1], c = 100, l = 65),
  chr_levels
)

coerce_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  map <- c(X = "23", x = "23", Y = "24", y = "24", MT = "25", Mt = "25", mt = "25", M = "25")
  x <- ifelse(x %in% names(map), map[x], x)
  suppressWarnings(as.integer(x))
}

prepare_gwas <- function(dt) {
  setDT(dt)
  if (!"POS" %in% names(dt)) setnames(dt, "GENPOS", "POS", skip_absent = TRUE)
  if (!"REF" %in% names(dt)) setnames(dt, "ALLELE0", "REF", skip_absent = TRUE)
  if (!"ALT" %in% names(dt)) setnames(dt, "ALLELE1", "ALT", skip_absent = TRUE)
  if (!"AF"  %in% names(dt)) setnames(dt, "A1FREQ", "AF", skip_absent = TRUE)
  if (!"CHR" %in% names(dt)) setnames(dt, "CHROM", "CHR", skip_absent = TRUE)

  dt[, CHR := coerce_chr(CHR)]
  if (autosomes_only) dt <- dt[CHR >= 1 & CHR <= 22]
  # ensure numeric to avoid int overflow later
  suppressWarnings({
    dt[, POS := as.numeric(POS)]
    dt[, P := as.numeric(P)]
  })
  dt <- dt[is.finite(P) & P > 0 & is.finite(POS)]
  dt[, `-log10P` := -log10(P)]

  # Build spans ensuring ALL chr 1..22 appear on axis (no skips)
  chr_span <- dt[, .(chr_len = max(POS, na.rm = TRUE)), by = CHR]
  setnames(chr_span, "CHR", "CHR_num")
  full_span <- data.table(CHR_num = 1:22)[chr_span, on = "CHR_num"]
  full_span[is.na(chr_len) | !is.finite(chr_len), chr_len := 1]              # keep tiny block for empty chr
  setorder(full_span, CHR_num)
  full_span[, cum_start := c(0, cumsum(as.numeric(chr_len))[-.N])]

  # Join back to data
  dt <- dt[full_span, on = .(CHR = CHR_num)]
  dt[, POS_cum := as.numeric(POS) + as.numeric(cum_start)]
  dt[, CHR_f := factor(as.character(CHR), levels = chr_levels)]

  # Axis ticks (centers)
  axis_dt <- full_span[, .(CHR = CHR_num, center = as.numeric(cum_start) + as.numeric(chr_len)/2)]
  attr(dt, "axis_dt") <- axis_dt
  dt
}

read_trait_gwas <- function(path_prefix, trait) {
  file_path <- file.path(path_prefix, paste0(trait, "_38_37.txt"))
  data.table::fread(file_path,
                    select = c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "P", "A1FREQ", "BETA", "SE"))
}

read_all_loci <- function(pathfile) {
  sep <- if (endsWith(tolower(pathfile), ".csv")) "," else "\t"
  al <- data.table::fread(pathfile, sep = sep)
  stopifnot("Phenotype" %in% names(al))
  if (!"CHROM" %in% names(al) && "CHR" %in% names(al)) setnames(al, "CHR", "CHROM")
  if (!"GENPOS" %in% names(al) && "POS" %in% names(al)) setnames(al, "POS", "GENPOS")
  al[, CHROM := coerce_chr(CHROM)]
  suppressWarnings(al[, `:=`(GENPOS = as.numeric(GENPOS))])
  if ("SNP" %in% names(al)) {
    parts <- tstrsplit(al$SNP, ":", fixed = TRUE)
    if (length(parts) >= 2) {
      al[, CHR_SNP := coerce_chr(parts[[1]])]
      suppressWarnings(al[, POS_SNP := as.numeric(parts[[2]])])
      if (length(parts) >= 4) {
        al[, REF_SNP := parts[[3]]]
        al[, ALT_SNP := parts[[4]]]
      }
    }
  }
  al
}

# pick lead SNPs separated by >= dist_kb within each chr, keep most significant per region
select_lead_labels <- function(dt, thr = genomewide, dist_kb = 500) {
  hits <- copy(dt)[P <= thr][order(P)]
  if (nrow(hits) == 0) return(hits[0])
  keep <- logical(nrow(hits))
  last_pos <- rep(-Inf, 22)
  for (i in seq_len(nrow(hits))) {
    cidx <- hits$CHR[i]
    if (is.na(cidx) || cidx < 1 || cidx > 22) next
    if (hits$POS[i] - last_pos[cidx] >= dist_kb * 1000) {
      keep[i] <- TRUE
      last_pos[cidx] <- hits$POS[i]
    }
  }
  hits[keep]
}

label_candidates_from_all_loci <- function(dt, trait, al, thr = genomewide, dist_kb = label_dist_kb) {
  lead <- select_lead_labels(dt, thr = thr, dist_kb = dist_kb)
  if (nrow(lead) == 0) return(lead)

  al_trait <- al[Phenotype == trait]
  if (nrow(al_trait) == 0) return(lead)

  # 1) exact CHR/GENPOS
  m1 <- merge(
    lead[, .(CHR, POS, P, POS_cum)],
    al_trait[, .(CHR = CHROM, POS = GENPOS, Locus_name)],
    by = c("CHR", "POS"),
    all.x = TRUE
  )

  # 2) fallback by parsed SNP CHR:POS using a clean table and then coalesce
  alt_tbl <- al_trait[!is.na(CHR_SNP) & !is.na(POS_SNP),
                      .(CHR = CHR_SNP, POS = POS_SNP, Locus_name_snp = Locus_name)]
  if (nrow(alt_tbl)) {
    m1 <- merge(m1, alt_tbl, by = c("CHR", "POS"), all.x = TRUE)
    m1[, Locus_name := fifelse(is.na(Locus_name), Locus_name_snp, Locus_name)]
    m1[, Locus_name_snp := NULL]
  }

  lab <- unique(m1[!is.na(Locus_name),
                   .(CHR, POS, POS_cum, P, label = Locus_name)])
  if (!nrow(lab)) return(lab[0])
  lab[, `-log10P` := -log10(P)]
  lab[, CHR_f := factor(as.character(CHR), levels = chr_levels)]
  lab
}

base_manhattan_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
}

plot_single_trait <- function(dt, trait, labels_dt = NULL, out_png, gw = genomewide, sug = suggestive) {
  axis_dt <- attr(dt, "axis_dt")
  dt[, CHR_f := factor(as.character(CHR), levels = chr_levels)]

  p <- ggplot(dt, aes(x = POS_cum, y = `-log10P`, color = CHR_f)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = chr_palette, drop = FALSE) +
    geom_hline(yintercept = -log10(gw), linetype = "dashed") +
    geom_hline(yintercept = -log10(sug), linetype = "dotted") +
    scale_x_continuous(
      breaks = axis_dt$center,
      labels = axis_dt$CHR,
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    labs(x = "Chromosome (1–22)", y = expression(-log[10](P)), title = trait) +
    base_manhattan_theme()

  if (!is.null(labels_dt) && nrow(labels_dt)) {
    p <- p +
      geom_point(data = labels_dt, aes(x = POS_cum, y = `-log10P`, color = CHR_f), size = 1.8) +
      ggrepel::geom_label_repel(
        data = labels_dt,
        aes(x = POS_cum, y = `-log10P`, label = label),
        size = 2.8, label.size = 0.15, segment.size = 0.2, max.overlaps = 100,
        min.segment.length = 0, box.padding = 0.25, point.padding = 0.2
      )
  }

  png(out_png, width = save_width_px, height = save_height_px, res = save_res)
  print(p)
  dev.off()
}

# ----------------------------- Execute -----------------------------
al <- read_all_loci(all_loci)

for (trait in traits) {
  gwas <- read_trait_gwas(gwas_path, trait) |> prepare_gwas()
  labs <- label_candidates_from_all_loci(gwas, trait, al, thr = genomewide, dist_kb = label_dist_kb)
  out <- file.path(plot_outpath, paste0("Manhattan_", trait, "_annotated.png"))
  plot_single_trait(gwas, trait, labels_dt = labs, out_png = out)
}
