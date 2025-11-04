library(tidyverse)
library(magrittr)
library(httr)
library(jsonlite)
library(xml2)
library(data.table)

# Source shared config (already in runroot from Nextflow)
source("config_R.R")

# CADD score retrieval via API
loci <- fread(var_file)

loci_lead <- loci %>% filter(SNP %in% unique(loci$lead_snp))

#loci_filtered <- loci %>% filter(credible_set_99 == 'Yes' | R2 >= 0.8)
loci_filtered <- loci
setnames(loci_filtered, c('ALLELE0', 'ALLELE1', 'GENPOS'), c('REF', 'ALT', 'POS'))
variants <- loci_filtered %>% select(SNP, ALT, REF, CHROM, POS) %>% mutate(POS = paste(CHROM, POS, sep = ':'))

variants_distinct <- variants %>% distinct(SNP, ALT, REF, .keep_all = TRUE)

url <- 'http://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/'

snp_to_query <- split(variants_distinct$POS, ceiling(seq_along(variants_distinct$POS) / 500))

total_iterations <- length(snp_to_query)  # Total number of iterations
print('Number of iterations:')
print(total_iterations)
iter <- 1
res_all <- list()

for (snp in snp_to_query) {
  for (i in snp) {
    r <- GET(paste(url, i, sep = ""), content_type("application/json"))
    res_all[[length(res_all) + 1]] <- fromJSON(toJSON(content(r)))
  }
  print(glue::glue("Iteration {iter} of {total_iterations} finished."))
  iter <- iter + 1
}


res_all_df <- res_all %>% 
  keep(~ is.data.frame(.x)) %>% 
  map(~ mutate_all(.x, unlist)) %>% 
  bind_rows(.id = 'SNP')

res_all_df <- res_all_df %>%
  mutate(SNP = paste(Chrom, Pos, Ref, Alt, sep = ":"))
  
res_all_df_output <- dplyr::select(res_all_df, -Alt, Chrom, Ref, Pos)
setnames(res_all_df_output, c('PHRED', 'RawScore'), c('CADD_PHRED', 'CADD_RawScore'))
res_all_df_output <- res_all_df_output %>%
  group_by(SNP) %>%
  slice_max(order_by = CADD_RawScore, n = 1, with_ties = FALSE) %>%
  ungroup()
fwrite(res_all_df_output, paste0(var_path, '/all_traits_in_ld_CADD.csv'))