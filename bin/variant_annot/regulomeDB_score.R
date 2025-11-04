library(data.table)
library(tidyverse)

source("config_R.R")

gwas <- fread(var_file)
print('Reading in regulomeDB data...')
regulomedb <- fread(paste0(databases, 'regulomedb/ENCFF250UJY_Dec23.tsv'))
regulomedb[, chrom := as.character(chrom)]

gwas$chrom <- paste0('chr', gwas$CHROM)
gwas$start <- gwas$GENPOS
gwas$end <- gwas$GENPOS
setkey(gwas, chrom, start, end)
setkey(regulomedb, chrom, start, end)

print('Finding overlapping SNPs...')
regulome_gwas <- foverlaps(gwas, regulomedb, nomatch = 0)

print('Number of overlapping SNPs:')
print(nrow(regulome_gwas))
regulome_gwas <- dplyr::select(regulome_gwas, SNP, ranking, probability_score)
regulome_gwas <- dplyr::rename(regulome_gwas, 'RegulomeDB_ranking' = ranking, 'RegulomeDB_probability_score' = probability_score)

regulome_gwas <- regulome_gwas[order(-RegulomeDB_probability_score), ]
regulome_gwas <- regulome_gwas[!duplicated(regulome_gwas$SNP), ]

gwas_out <- merge(gwas, regulome_gwas, by = 'SNP', all.x = TRUE)
gwas_out$chrom <- NULL
gwas_out$start <- NULL
gwas_out$end <- NULL
gwas_out <- dplyr::select(gwas_out, SNP, RegulomeDB_ranking, RegulomeDB_probability_score)

gwas_out <- gwas_out %>%
  group_by(SNP) %>%
  slice_max(order_by = RegulomeDB_probability_score, n = 1, with_ties = FALSE) %>%
  ungroup()

fwrite(gwas_out, paste0(var_path, '/all_traits_in_ld_regulomedb.txt'), sep = '\t')
print("Saved output data:")
print(head(gwas_out))