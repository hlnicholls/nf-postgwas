#!/usr/bin/env Rscript

# Enrichr gene annotation (robust, single-threaded)
suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(data.table)
  library(enrichR)
})

try({
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    R_DEFAULT_INTERNET_TIMEOUT = "120"
  )
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
  options(mc.cores = 1)
  # Silence harmless systemd warning in containers
  try(system("timedatectl 2>/dev/null", intern = TRUE), silent = TRUE)
}, silent = TRUE)

# ---- Load NF-shimmed config (defines var_file, enrichment_path, etc.) ----
source("config_R.R")

# ---- Read & clean input genes (expects Gene_Symbol in var_file) ----
genes_all <- fread(var_file)
genes_all <- genes_all %>%
  dplyr::rename(Gene = Gene_Symbol) %>%
  dplyr::select(Gene) %>%
  dplyr::filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene, .keep_all = TRUE)

if (!nrow(genes_all)) {
  stop("No genes available after cleaning 'Gene_Symbol' from: ", var_file)
}

# ---- Discover available Enrichr libraries and intersect with desired set ----
db_list <- tryCatch(listEnrichrDbs(), error = function(e) NULL)
available <- if (!is.null(db_list) && "libraryName" %in% names(db_list)) {
  db_list$libraryName
} else character(0)

# Your requested/legacy set (kept intact)
desired <- c(
  "Chromosome_Location", "Human_Phenotype_Ontology",
  "OMIM_Expanded", "Disease_Signatures_from_GEO_down_2014",
  "Disease_Perturbations_from_GEO_down",
  "Disease_Perturbations_from_GEO_up", "Drug_Perturbations_from_GEO_down",
  "Drug_Perturbations_from_GEO_up", "Gene_Perturbations_from_GEO_up",
  "Gene_Perturbations_from_GEO_down", "Reactome_2022", "BioCarta_2016",
  "Panther_2016", "Jensen_TISSUES", "Jensen_COMPARTMENTS", "Jensen_DISEASES",
  "DSigDB", "HDSigDB_Human_2021", "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "dbGaP",
  "GWAS_Catalog_2023", "WikiPathway_2023_Human", "KEGG_2021_Human",
  "InterPro_Domains_2019", "MGI_Mammalian_Phenotype_Level_4_2021",
  "UK_Biobank_GWAS_v1", "BioPlanet_2019", "ClinVar_2019", "PheWeb_2019",
  "MGI_Mammalian_Phenotype_Level_3", "DisGeNET", "Elsevier_Pathway_Collection",
  "Proteomics_Drug_Atlas_2023", "The_Kinase_Library_2023",
  "GTEx_Tissues_V8_2023", "MAGNET_2023"
)

dbs <- if (length(available)) {
  intersect(desired, available)
} else {
  # Fallback if library list cannot be fetched
  c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022")
}

if (!length(dbs)) {
  stop("No overlap between desired Enrichr libraries and what is available on the server.")
}

# ---- Safe wrapper around enrichr() returning standardized tibble rows per DB ----
safe_enrichr_one <- function(gene, dbs) {
  res <- tryCatch(enrichr(gene, dbs), error = function(e) NULL)
  out <- vector("list", length(dbs)); names(out) <- dbs
  if (is.null(res)) return(out)  # all NULLs if failed

  for (d in dbs) {
    if (!is.null(res[[d]]) && nrow(res[[d]]) > 0) {
      db_table <- as_tibble(res[[d]])
      standardized_table <- db_table %>%
        mutate(
          Term            = if ("Term" %in% names(db_table)) as.character(Term) else NA_character_,
          Overlap         = if ("Overlap" %in% names(db_table)) as.character(Overlap) else NA_character_,
          P.value         = if ("P.value" %in% names(db_table)) as.numeric(P.value) else NA_real_,
          Adjusted.P.value= if ("Adjusted.P.value" %in% names(db_table)) as.numeric(Adjusted.P.value) else NA_real_,
          Odds.Ratio      = if ("Odds.Ratio" %in% names(db_table)) as.numeric(Odds.Ratio) else NA_real_,
          Combined.Score  = if ("Combined.Score" %in% names(db_table)) as.numeric(Combined.Score) else NA_real_,
          Genes           = if ("Genes" %in% names(db_table)) as.character(Genes) else NA_character_,
          Gene            = gene,
          Database        = d
        )
      out[[d]] <- standardized_table
    } else {
      out[[d]] <- NULL
    }
  }
  out
}

# ---- Run sequentially (no mclapply) to avoid OOM and OOB indexing ----
results_list <- vector("list", nrow(genes_all))
names(results_list) <- genes_all$Gene

for (i in seq_len(nrow(genes_all))) {
  g <- genes_all$Gene[i]
  message(sprintf("[%d/%d] %s", i, nrow(genes_all), g))
  results_list[[i]] <- safe_enrichr_one(g, dbs)
  if ((i %% 50) == 0) gc(FALSE)
}

# ---- Collate all non-null DB tables into one long table ----
rows <- list()
idx <- 1L
for (i in seq_along(results_list)) {
  gene_map <- results_list[[i]]
  if (is.null(gene_map)) next
  gene_frames <- Filter(Negate(is.null), gene_map)
  if (!length(gene_frames)) next
  rows[[idx]] <- bind_rows(gene_frames)
  idx <- idx + 1L
}

if (!length(rows)) {
  warning("Enrichr returned no results for any gene/database; writing empty table.")
  # Keep filename pattern the same as your original pipeline
  out_file <- file.path(enrichment_path, paste0("All_genes_annotated_with_EnrichR_", Sys.Date(), ".csv"))
  dir.create(enrichment_path, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(`Gene symbol` = character()), out_file)
  message("Wrote (empty): ", out_file)
  quit(save = "no", status = 0)
}

results_df <- bind_rows(rows)

# ---- Build the wide, per-gene table exactly like your original ----
filtered_table <- results_df %>%
  select(Gene, Database, Term) %>%
  mutate(Term = ifelse(is.na(Term) | Term == "", NA_character_, Term)) %>%
  group_by(Gene, Database) %>%
  summarise(Term = paste(na.omit(unique(Term)), collapse = ";"), .groups = "drop") %>%
  pivot_wider(
    names_from = Database,
    values_from = Term,
    names_sep = "_",
    values_fill = NA_character_
  ) %>%
  rename(`Gene symbol` = Gene)

# ---- Write output where NF expects it ----
dir.create(enrichment_path, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(enrichment_path, paste0("All_genes_annotated_with_EnrichR_", Sys.Date(), ".csv"))
fwrite(filtered_table, out_file)
message("Wrote: ", out_file)
