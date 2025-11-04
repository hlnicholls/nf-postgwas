#!/usr/bin/env Rscript

# Identify top PoPS genes per locus

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(topr)
  library(purrr)
  library(here)
})


source("config_R.R")

pops <- fread(paste0(pops_output_path,'/pops_results_all_features_cleaned.csv'))
colnames(pops)[1] <- 'Gene_Symbol'
pops <- filter(pops, Gene_Symbol != "")

for (phenotype in traits) {

  pops_columns <- grep(paste0(phenotype, "_pops_"), names(pops), value = TRUE)
  pops <- pops %>%
    rowwise() %>%
    mutate(!!paste0(phenotype, "_Highest_PoPS") := max(c_across(all_of(pops_columns)), na.rm = TRUE)) %>%
    ungroup()
}


gwas_in_ld <- fread(var_file,
                    select=c('Locus_name','lead_snp', 'Phenotype','Gene_Symbol'))


df <- merge(pops, gwas_in_ld, by='Gene_Symbol', all.x=TRUE, allow.cartesian=TRUE)

locus_groups <- fread(locus_blocks_output, select=c('Lead_SNP', 'Locus_number'))
colnames(locus_groups)[1] <- 'lead_snp'
n_unique_snps <- n_distinct(locus_groups$lead_snp)
df <- merge(df, locus_groups, by='lead_snp', all.x=TRUE)
n_unique_snps <- n_distinct(df$lead_snp)
df <- filter(df, !is.na(df$Gene_Symbol) & Gene_Symbol != '')
df[df == -Inf] <- NA

pops_columns <- grep("_Highest_PoPS", names(df), value = TRUE)

process_pops_column <- function(column_name) {
  df %>%
    group_by(Locus_number) %>%
    slice(which.max(get(column_name))) %>%
    select(lead_snp, Gene_Symbol, all_of(column_name))
}

results <- lapply(pops_columns, process_pops_column)

df2 <- df %>%
  mutate(Max_pops_value = pmap_dbl(select(., all_of(pops_columns)), ~ {
    values <- c(...)
    if (all(is.na(values))) {
      return(NA_real_)
    } else {
      return(max(values, na.rm = TRUE))
    }
  }))


selected_genes <- df2 %>%
  group_by(Locus_number) %>%
  slice_max(Max_pops_value, n = 1, with_ties = FALSE) %>%
  ungroup()

selected_genes <- selected_genes %>% rename(Nearest_Gene_10kb = Gene_Symbol)

final_df <- selected_genes %>%
  dplyr::select(Locus_number, Locus_name, Nearest_Gene_10kb, all_of(pops_columns))

df$Top_PoPS_Score_per_Locus <- ifelse(df$Gene_Symbol %in% final_df$Nearest_Gene_10kb, 'Yes', 'No')
df <- df %>% rename(Nearest_Gene_10kb = Gene_Symbol)
output <- dplyr::select(df, Nearest_Gene_10kb, Top_PoPS_Score_per_Locus, Locus_name, Locus_number)
genes <- as.data.frame(unique(selected_genes$Nearest_Gene_10kb))
colnames(genes)[1] <- 'Gene'

fwrite(genes, paste0(pops_output_path,'/pops_top_genes_per_locus.txt'))

fwrite(output, paste0(pops_output_path,'/gwas_all_loci_top_pops_genes.txt'))

top_pops <- filter(df, Top_PoPS_Score_per_Locus == 'Yes')
top_pops <- dplyr::select(top_pops, Nearest_Gene_10kb, Top_PoPS_Score_per_Locus)
top_pops <- unique(top_pops)
cat(sprintf("Top PoPS genes identified: %d\n", nrow(top_pops)))

cat("PoPS per loci analysis completed successfully\n")