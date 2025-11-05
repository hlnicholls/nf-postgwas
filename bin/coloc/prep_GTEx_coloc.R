# eQTL colocalisation analysis pipeline

# Source the runroot config file explicitly from the current working directory.
# The Nextflow process changes into the runroot before running this script, so
# sourcing './config_R.R' is the most robust approach. Fail early with a clear
# error if the file is missing to avoid vague "object not found" errors later.
cfg <- file.path(getwd(), "config_R.R")
if (!file.exists(cfg)) {
  stop(sprintf("Missing runroot config: %s\nEnsure CREATE_CONFIG_SHIMS ran and produced config_R.R in the runroot directory.", cfg))
}
source(cfg)
library(data.table)

output_dir <- paste0(coloc_path,"/eQTL")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (phenotype in traits){

  data.ld38 = fread(var_file)
  
  setwd(eqtl_in)
  
  data.ld38 = subset(data.ld38, R2>=0.8) 
  data.ld38$ID = seq_len(nrow(data.ld38))
  
  data.ld38$variant_id = paste("chr",data.ld38$CHROM,"_",data.ld38$GENPOS,"_",data.ld38$ALLELE1,"_",data.ld38$ALLELE0,"_b38", sep="")
  data.ld38$variant_id2 = paste("chr",data.ld38$CHROM,"_",data.ld38$GENPOS,"_",data.ld38$ALLELE0,"_",data.ld38$ALLELE1,"_b38", sep="")
  
  cs.egene.all <- data.frame()
  
  for (tissue in  tissues){
    
    # read significant egenes from GTEX 
    egene = read.table(paste(gtex_path, tissue,".v8.egenes.txt.gz", sep=""),h=T,stringsAsFactors = F, sep="\t")
    # Note that the *.egenes.txt.gz files contain data for all genes tested; to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
    # selected significant genes (as instructed by GTEx) based on Storey q-value derived from pval_beta
    egene = subset(egene, qval<=0.05)
    
    cs.egene = merge(egene, data.ld38, by = "variant_id")
    setnames(data.ld38, old = c("variant_id", "variant_id2"), new = c("variant_id2", "variant_id"))
    cs.egene2 = merge(egene, data.ld38, by = "variant_id")
    setnames(data.ld38, old = c("variant_id2", "variant_id"), new = c("variant_id", "variant_id2"))
    cs.egene = rbind(cs.egene,cs.egene2)
    
    if(nrow(cs.egene)>0){
      cs.egene$tissue = tissue
      cs.egene.all = rbind(cs.egene.all,cs.egene)
    }
    
  }
  
  # these are the candidate that we will verify using colocalisation analysis
  write.table(cs.egene.all,paste0(phenotype,"_eQTL_candidates.txt"), col.names = T, sep="\t")
}

