#!/usr/bin/env Rscript

# IMPC enrichment analysis for prioritised genes
# Uses disease_terms / disease_regex provided by CREATE_CONFIG_SHIMS (config_R.R)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GeneOverlap)
  library(splitstackshape)
})

# ---------- Utilities ----------
safe_exists <- function(sym, env = .GlobalEnv) exists(sym, envir = env, inherits = FALSE)
safe_get    <- function(sym, default = NULL, env = .GlobalEnv) if (safe_exists(sym, env)) get(sym, envir = env, inherits = FALSE) else default
nz_or_blank <- function(x) ifelse(is.na(x), "", x)

write_empty_outputs <- function(out_dir, why = "No data") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Ranked ORA table
  ranked_path <- file.path(out_dir, "IMPC_genes_ranked_enriched.csv")
  fwrite(
    data.frame(Phenotype = character(), tested = integer(), pval = numeric(), pval_adjusted = numeric()),
    ranked_path
  )

  # Collapsed results table
  collapsed_path <- file.path(out_dir, "IMPC_results_table.csv")
  fwrite(
    data.frame(Phenotype = character(), pval = numeric(), pval_adjusted = numeric(),
               `MGI Gene Id` = character(), Gene = character()),
    collapsed_path
  )

  # Placeholder plot
  plot_path <- file.path(out_dir, "Genes_ranked_IMPC_enrichment.png")
  png(plot_path, width = 14, height = 8, units = "in", res = 300)
  plot.new()
  title(main = paste("IMPC enrichment:", why))
  dev.off()
}

# ---------- Config ----------
cfg_loaded <- FALSE
try({
  source("config_R.R", local = TRUE)
  cfg_loaded <- TRUE
}, silent = TRUE)

# Paths (from config when available)
enrichment_output_path <- tryCatch(get("enrichment_output_path"), error = function(e) getwd())
databases              <- tryCatch(get("databases"),              error = function(e) NULL)

# Prioritised genes path (from config variable 'prioritised_genes' which is a file path)
prioritised_genes_path <- tryCatch(get("prioritised_genes"), error = function(e) NULL)
if (is.null(prioritised_genes_path) || is.na(prioritised_genes_path) || prioritised_genes_path == "") {
  # Fallback to a local file if the config shim wasn't created (keeps your wrapper happy)
  prioritised_genes_path <- "Prioritised_genes.csv"
}

dir.create(enrichment_output_path, showWarnings = FALSE, recursive = TRUE)

# ---------- Disease terms / regex with safe default ----------
default_disease_terms <- c(
  "cardio","cardiac","atrial","myocard","arrhyt","vascular","vessel","vasculature",
  "heart","hypertroph","dilated","ventric","hypertension","blood pressure","diabetes",
  "obesity","hypotension","qrs","interval","segment","aorta","cholesterol"
)
disease_terms  <- safe_get("disease_terms", character())
disease_regex  <- safe_get("disease_regex", NULL)

if (is.null(disease_regex) || isTRUE(nchar(nz_or_blank(disease_regex)) == 0)) {
  if (length(disease_terms)) {
    disease_regex <- paste(disease_terms, collapse = "|")
  } else {
    disease_regex <- paste(default_disease_terms, collapse = "|")
    warning("IMPC_enrichment: No disease terms found in config_R.R — using default regex.")
  }
}
message("IMPC_enrichment: using disease regex: ", disease_regex)

# ---------- Load inputs (prioritised genes) ----------
if (!file.exists(prioritised_genes_path) || isTRUE(file.info(prioritised_genes_path)$size == 0)) {
  message("IMPC_enrichment: Prioritised genes CSV missing/empty: ", prioritised_genes_path)
  write_empty_outputs(enrichment_output_path, "No prioritised genes")
  quit(save = "no", status = 0)
}

pri <- tryCatch(fread(prioritised_genes_path, na.strings = c("", "NA")),
                error = function(e) data.table())
if (nrow(pri) == 0) {
  message("IMPC_enrichment: Prioritised genes CSV has 0 rows: ", prioritised_genes_path)
  write_empty_outputs(enrichment_output_path, "No prioritised genes")
  quit(save = "no", status = 0)
}

# Ensure essential columns exist (fill if missing)
if (!("Nearest_Gene_10kb" %in% names(pri))) pri$Nearest_Gene_10kb <- NA_character_
if (!("Gene_Prioritisation_Score" %in% names(pri))) pri$Gene_Prioritisation_Score <- NA_real_

top_pops <- pri %>%
  select(Nearest_Gene_10kb, Gene_Prioritisation_Score) %>%
  filter(!is.na(Gene_Prioritisation_Score), Gene_Prioritisation_Score != 0) %>%
  distinct(Nearest_Gene_10kb) %>%
  rename(Gene = Nearest_Gene_10kb) %>%
  mutate(Gene = toupper(Gene)) %>%
  filter(!is.na(Gene) & Gene != "")

prior_genes <- unique(top_pops$Gene)
if (length(prior_genes) == 0) {
  message("IMPC_enrichment: No prioritised genes with non-zero scores.")
  write_empty_outputs(enrichment_output_path, "No prioritised genes with non-zero scores")
  quit(save = "no", status = 0)
}

# ---------- Load IMPC table ----------
impc_path <- if (!is.null(databases)) file.path(databases, "impc", "phenotypeHitsPerGene_20250416.csv") else NA_character_
if (is.na(impc_path) || !file.exists(impc_path) || isTRUE(file.info(impc_path)$size == 0)) {
  message("IMPC_enrichment: IMPC file missing/empty: ", impc_path)
  write_empty_outputs(enrichment_output_path, "IMPC source missing/empty")
  quit(save = "no", status = 0)
}

mouse <- tryCatch(fread(impc_path, na.strings = c("", "NA")),
                  error = function(e) data.table())
if (nrow(mouse) == 0) {
  message("IMPC_enrichment: IMPC table has 0 rows: ", impc_path)
  write_empty_outputs(enrichment_output_path, "Empty IMPC table")
  quit(save = "no", status = 0)
}

# Normalise gene symbol column name
if ("Gene Symbol" %in% names(mouse)) setnames(mouse, "Gene Symbol", "Gene")
if (!("Gene" %in% names(mouse))) mouse$Gene <- NA_character_
mouse <- mouse %>% mutate(Gene = toupper(Gene))

# Ensure phrase columns exist (stay compatible with previous expectations)
if (!("# Phenotype Hits" %in% names(mouse))) mouse[["# Phenotype Hits"]] <- NA_integer_
if (!("Phenotype Hits" %in% names(mouse)))  mouse[["Phenotype Hits"]]     <- NA_character_

# ---------- Build full and disease-related subsets ----------
mouse_all <- top_pops %>%
  left_join(mouse, by = "Gene") %>%
  distinct()

# Disease-related rows: any disease term in string
mouse_disease <- mouse_all %>%
  filter(!is.na(`Phenotype Hits`) & grepl(disease_regex, `Phenotype Hits`, ignore.case = TRUE))

# ---------- Split phenotypes into long form ----------
to_long <- function(df) {
  if (nrow(df) == 0) {
    # return an empty table with expected columns
    return(data.frame(Gene = character(), Phenotype = character(), stringsAsFactors = FALSE))
  }
  if (!("Phenotype Hits" %in% names(df))) {
    return(data.frame(Gene = df$Gene %||% character(0), Phenotype = character(0), stringsAsFactors = FALSE))
  }
  out <- tryCatch(
    cSplit(df, "Phenotype Hits", sep = "::", direction = "long"),
    error = function(e) {
      # fallback: treat entire string as one phenotype
      data.frame(Gene = df$Gene, `Phenotype Hits` = df[["Phenotype Hits"]], stringsAsFactors = FALSE)
    }
  )
  setnames(out, "Phenotype Hits", "Phenotype")
  out
}
`%||%` <- function(a,b) if (!is.null(a)) a else b

mouse_all_long     <- to_long(mouse_all)     %>% select(Gene, Phenotype)
mouse_disease_long <- to_long(mouse_disease) %>% select(Gene, Phenotype)

# Remove NA/blank phenotype labels
mouse_all_long     <- mouse_all_long     %>% filter(!is.na(Phenotype), Phenotype != "", !is.na(Gene), Gene != "")
mouse_disease_long <- mouse_disease_long %>% filter(!is.na(Phenotype), Phenotype != "", !is.na(Gene), Gene != "")

# Universe: all genes with any IMPC phenotype listed
impc_universe <- unique(mouse$Gene[!is.na(mouse$Gene) & mouse$Gene != ""])
genome_size   <- max(length(impc_universe), 20000)  # keep a large stable universe if needed

# ---------- ORA: per-phenotype overlap tests ----------
run_ora <- function(prior_genes, impc_long, universe_size) {
  if (nrow(impc_long) == 0 || length(prior_genes) == 0) {
    return(data.frame(Phenotype = character(), tested = integer(), pval = numeric(), pval_adjusted = numeric()))
  }
  split_impc <- split(impc_long$Gene, impc_long$Phenotype)
  tests <- lapply(split_impc, function(genes_in_pheno) {
    genes_in_pheno <- unique(genes_in_pheno[!is.na(genes_in_pheno) & genes_in_pheno != ""])
    go.obj <- newGeneOverlap(prior_genes, genes_in_pheno, genome.size = universe_size)
    res    <- testGeneOverlap(go.obj)
    list(tested = getTested(res), pval = getPval(res))
  })
  res <- data.frame(
    Phenotype     = names(tests),
    tested        = vapply(tests, `[[`, integer(1), "tested"),
    pval          = vapply(tests, `[[`, numeric(1),  "pval"),
    stringsAsFactors = FALSE
  )
  res$pval_adjusted <- p.adjust(res$pval, "fdr")
  res[order(res$pval_adjusted, decreasing = FALSE), ]
}

ora_all     <- run_ora(prior_genes, mouse_all_long, genome_size)
ora_disease <- run_ora(prior_genes, mouse_disease_long, genome_size) # kept for reference

# ---------- Save main table (ranked by adjusted p) ----------
ranked_path <- file.path(enrichment_output_path, "IMPC_genes_ranked_enriched.csv")
if (nrow(ora_all) == 0) {
  fwrite(
    data.frame(Phenotype = character(), tested = integer(), pval = numeric(), pval_adjusted = numeric()),
    ranked_path
  )
} else {
  fwrite(ora_all, ranked_path)
}

# ---------- Summaries for plotting (top 20 by count among prioritised genes) ----------
mouse_hits_for_prior <- mouse_all_long %>%
  filter(Gene %in% prior_genes) %>%
  count(Phenotype, name = "Count") %>%
  arrange(desc(Count))

top20 <- head(mouse_hits_for_prior, 20)

plot_df <- top20 %>%
  left_join(ora_all %>% select(Phenotype, pval_adjusted), by = "Phenotype") %>%
  mutate(
    pval_adjusted = signif(pval_adjusted, 5),
    Phenotype = factor(Phenotype, levels = Phenotype)
  )

plot_path <- file.path(enrichment_output_path, "Genes_ranked_IMPC_enrichment.png")
png(plot_path, width = 14, height = 8, units = "in", res = 300)
if (nrow(plot_df) == 0) {
  plot.new()
  title("Top 20 Most Frequent IMPC Mouse Phenotypes (no data)")
} else {
  # Compute a label position slightly inside the bar end (works with coord_flip)
  max_count <- max(plot_df$Count, na.rm = TRUE)
  plot_df <- plot_df %>%
    mutate(
      label_y = ifelse(is.na(Count), NA_real_, pmax(0, Count - 0.02 * max_count)),
      label_text = ifelse(is.na(pval_adjusted), "", as.character(pval_adjusted))
    )

  ggplot(plot_df, aes(x = Phenotype, y = Count)) +
    geom_bar(stat = "identity", fill = "#5e2a84") +   # deep purple (kept)
    coord_flip() +
    ggtitle("Top 20 Most Frequent IMPC Mouse Phenotypes") +
    ylab("Gene Count") +
    xlab("Phenotype") +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      axis.title   = element_text(size = 14),
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    # Adjusted p-values shown IN WHITE on top (near the end) of each bar
    geom_text(aes(y = label_y, label = label_text), color = "white", size = 3.5, na.rm = TRUE) %>%
    print()
}
dev.off()

# ---------- Collapsed results (match previous shape) ----------
# If MGI Gene Id exists in source, join and collapse; otherwise keep empty string column
if ("MGI Gene Id" %in% names(mouse)) {
  mouse_all_long_ids <- mouse_all_long %>%
    left_join(mouse[, c("Gene", "MGI Gene Id")], by = "Gene") %>%
    distinct()
  collapsed <- ora_all %>%
    left_join(mouse_all_long_ids, by = "Phenotype") %>%
    group_by(Phenotype) %>%
    summarise(
      pval          = dplyr::first(pval),
      pval_adjusted = dplyr::first(pval_adjusted),
      `MGI Gene Id` = toString(unique(`MGI Gene Id`[!is.na(`MGI Gene Id`)])),
      Gene          = toString(unique(Gene)),
      .groups = "drop"
    )
} else {
  collapsed <- ora_all %>%
    left_join(mouse_all_long, by = "Phenotype") %>%
    group_by(Phenotype) %>%
    summarise(
      pval          = dplyr::first(pval),
      pval_adjusted = dplyr::first(pval_adjusted),
      `MGI Gene Id` = "",
      Gene          = toString(unique(Gene)),
      .groups = "drop"
    )
}

collapsed_path <- file.path(enrichment_output_path, "IMPC_results_table.csv")
if (nrow(collapsed) == 0) {
  fwrite(
    data.frame(Phenotype = character(), pval = numeric(), pval_adjusted = numeric(),
               `MGI Gene Id` = character(), Gene = character()),
    collapsed_path
  )
} else {
  fwrite(collapsed, collapsed_path)
}

message("IMPC enrichment analysis completed successfully")
