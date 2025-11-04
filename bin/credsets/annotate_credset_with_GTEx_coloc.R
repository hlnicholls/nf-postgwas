suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

source("config_R.R")

credsets <- fread(credidble_sets_path)
coloc_all <- data.frame()

for (trait in traits) {
  coloc_res <- fread(paste0(coloc_path, '/', trait, '_eQTL_COLOC.tsv'), select=c("SNP.x", "tissue", "PP3", "PP4"))
  coloc_res <- coloc_res %>% rename(SNP = SNP.x,
                                      GTEx_coloc_tissue = tissue,
                                      GTEx_coloc_PP3 = PP3,
                                      GTEx_coloc_PP4 = PP4,
                                      )
  coloc_res$Phenotype <- trait
  coloc_res <- distinct(coloc_res)

  collapse_coloc <- function(df) {
    df[, PP4 := as.numeric(GTEx_coloc_PP4)]
    df[, PP3 := as.numeric(GTEx_coloc_PP3)]
    df[, tissue_fmt := paste0(GTEx_coloc_tissue, " (PP3: ", PP3, ") (PP4: ", PP4, ")")]
    
    df[, .(
      GTEx_coloc_tissue = paste(tissue_fmt, collapse = ", "),
      has_GTEx_coloc_PP4_08 = ifelse(any(PP4 >= 80), "Yes", "No"),
      GTEx_coloc_PP4_08_Tissue = paste(GTEx_coloc_tissue[PP4 >= 80], collapse = ", ")
    ), by = .(SNP, Phenotype)]
  }

  coloc_res <- collapse_coloc(coloc_res)

  coloc_all <- rbind(coloc_all, coloc_res)
}


credsets <- merge(credsets, coloc_all, all.x = TRUE)
credsets$has_GTEx_coloc_PP4_08[is.na(credsets$has_GTEx_coloc_PP4_08)] <- "No"


fwrite(credsets, paste0(output_path, '/Credible_set_annotation','/Annotated_credsets_LD08_GTExColoc.csv'))

