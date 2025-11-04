# bin/credsets/GWAS_Catalog_Query.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(gwasrapidd)
  library(data.table)
})

source("config_R.R") 
# -------- helpers --------
write_empty_and_quit <- function(outdir, filename_prefix = "GWAS_catalog_query_GWAS_RAPIDD_") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # final expected columns (names after renaming)
  target_cols <- c(
    "Locus_name","Phenotype","Method","CHR","Lead_variant","SNP","R2",
    "BP","EA","NEA","BETA","SE","P",
    "GWAS_Catalog_DISEASE/TRAIT","GWAS_Catalog_EFO_TRAIT",
    "GWAS_Catalog_INITIAL SAMPLE SIZE","GWAS_Catalog_REPLICATION SAMPLE SIZE",
    "GWAS_Catalog_GxE","GWAS_Catalog_GxG","GWAS_Catalog_RISK ALLELE",
    "GWAS_Catalog_RISK Frequency","GWAS_Catalog_multiple_snp_haplotype",
    "GWAS_Catalog_BETA","GWAS_Catalog_OR","GWAS_Catalog_SE",
    "GWAS_Catalog_95% CI","GWAS_Catalog_P","GWAS_Catalog_BETA_unit",
    "GWAS_Catalog_BETA_direction","GWAS_Catalog_BETA_description","GWAS_Catalog_PMID"
  )

  empty_df <- as_tibble(setNames(rep(list(vector(mode = "character")), length(target_cols)), target_cols))
  outfile <- file.path(outdir, paste0(filename_prefix, Sys.Date(), ".csv"))
  data.table::fwrite(empty_df, outfile)
  quit(save = "no", status = 0)
}

ensure_cols <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    for (nm in missing) df[[nm]] <- NA
  }
  df
}

safe_bind_traits <- function(assoc_ids) {
  if (length(assoc_ids) == 0) {
    return(tibble(association_id = character(), trait = character(), reported_trait = character()))
  }
  tr_list <- purrr::map(assoc_ids, function(x) {
    tryCatch({
      tr <- get_traits(association_id = x)
      # ensure we always have a traits table
      tr@traits %>% mutate(association_id = x)
    }, error = function(e) {
      tibble(association_id = x, trait = NA_character_, reported_trait = NA_character_)
    })
  })
  bind_rows(tr_list)
}

# -------- load inputs --------
var_annots <- data.table::fread(credidble_sets_path)
lead_var   <- data.table::fread(all_loci)

# explode rsIDs (comma-separated) & clean
var_annots_distinct <- var_annots %>%
  distinct(rsid_1kg, .keep_all = TRUE) %>%
  separate_rows(rsid_1kg, sep = ",") %>%
  mutate(rsid_1kg = trimws(rsid_1kg)) %>%
  filter(rsid_1kg != "")

# If no rsIDs at all -> write empty and exit 0
if (nrow(var_annots_distinct) == 0) {
  write_empty_and_quit(pleiotropy_outpath)
}

# -------- query GWAS Catalog --------
ga <- tryCatch({
  get_associations(variant_id = var_annots_distinct$rsid_1kg, verbose = TRUE)
}, error = function(e) {
  # If the API barfs, treat as no results but do not fail the pipeline.
  NULL
})

if (is.null(ga)) {
  write_empty_and_quit(pleiotropy_outpath)
}

associations <- ga@associations
alleles      <- ga@risk_alleles

# If nothing came back, exit cleanly
if (nrow(associations) == 0) {
  write_empty_and_quit(pleiotropy_outpath)
}

# Ensure join keys exist
associations <- ensure_cols(associations, c("association_id"))
alleles      <- ensure_cols(alleles,      c("association_id"))

assoc_with_alleles <- associations %>% left_join(alleles, by = "association_id")

# Traits (guard per-association)
traits <- safe_bind_traits(assoc_with_alleles$association_id)

# make sure the join column is present even if traits is empty
traits <- ensure_cols(traits, c("association_id"))

assoc_with_traits <- assoc_with_alleles %>% left_join(traits, by = "association_id") %>% distinct(association_id, .keep_all = TRUE)

# If still nothing usable, exit cleanly
if (nrow(assoc_with_traits) == 0) {
  write_empty_and_quit(pleiotropy_outpath)
}

# Studies & pubs (guard everything)
assoc_ids <- assoc_with_traits$association_id

studies_association_map <- tryCatch(association_to_study(association_id = assoc_ids),
                                    error = function(e) NULL)
studies <- tryCatch(get_studies(association_id = assoc_ids),
                    error = function(e) NULL)

if (is.null(studies) || is.null(studies_association_map)) {
  study_ids <- tibble(study_id = character(), association_id = character())
  pubmed    <- tibble(study_id = character(), association_id = character(), pubmed_id = character())
} else {
  study_ids <- studies@studies %>%
    left_join(studies_association_map, by = "study_id") %>%
    ensure_cols(c("association_id"))

  pubmed <- studies@publications %>%
    left_join(studies_association_map, by = "study_id") %>%
    ensure_cols(c("association_id", "pubmed_id"))
}

assoc_with_studies <- assoc_with_traits %>%
  left_join(study_ids, by = "association_id") %>%
  left_join(pubmed,    by = "association_id") %>%
  distinct()

# Link back to our annotated credible set table
credible_overlap_tbl <- var_annots_distinct %>%
  inner_join(assoc_with_studies, by = c("rsid_1kg" = "variant_id"))

# If no overlap with our credible set rsIDs, exit empty
if (nrow(credible_overlap_tbl) == 0) {
  write_empty_and_quit(pleiotropy_outpath)
}

# Column selection/renaming (tolerant)
cols_to_choose <- c(
  "Locus_name","Phenotype","Method","CHROM","lead_snp","rsid_1kg","R2","GENPOS",
  "ALLELE1","ALLELE0","BETA","SE","P",
  "trait","reported_trait","initial_sample_size","replication_sample_size",
  "gxe","gxg","risk_allele","risk_frequency","multiple_snp_haplotype",
  "beta_number","or_per_copy_number","standard_error","range","pvalue",
  "beta_unit","beta_direction","beta_description","pubmed_id"
)
new_names <- c(
  "Locus_name","Phenotype","Method","CHR","Lead_variant","SNP","R2","BP",
  "EA","NEA","BETA","SE","P",
  "GWAS_Catalog_DISEASE/TRAIT","GWAS_Catalog_EFO_TRAIT",
  "GWAS_Catalog_INITIAL SAMPLE SIZE","GWAS_Catalog_REPLICATION SAMPLE SIZE",
  "GWAS_Catalog_GxE","GWAS_Catalog_GxG","GWAS_Catalog_RISK ALLELE",
  "GWAS_Catalog_RISK Frequency","GWAS_Catalog_multiple_snp_haplotype",
  "GWAS_Catalog_BETA","GWAS_Catalog_OR","GWAS_Catalog_SE","GWAS_Catalog_95% CI",
  "GWAS_Catalog_P","GWAS_Catalog_BETA_unit","GWAS_Catalog_BETA_direction",
  "GWAS_Catalog_BETA_description","GWAS_Catalog_PMID"
)

credible_overlap_tbl <- ensure_cols(credible_overlap_tbl, cols_to_choose)
credible_overlap_tbl <- credible_overlap_tbl %>% dplyr::select(all_of(cols_to_choose))
names(credible_overlap_tbl) <- new_names

credible_overlap_distinct <- credible_overlap_tbl %>% distinct(.keep_all = TRUE)

dir.create(pleiotropy_outpath, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(pleiotropy_outpath, paste0("GWAS_catalog_query_GWAS_RAPIDD_", Sys.Date(), ".csv"))
data.table::fwrite(credible_overlap_distinct, outfile)
