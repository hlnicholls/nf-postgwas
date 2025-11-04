library(data.table)

source("config_R.R")

gwas_data <- fread(var_file)
gwas_data <- gwas_data[order(gwas_data$CHROM, as.numeric(gwas_data$GENPOS)), ]
gwas_data$FORMAT <- NA
gwas_data$CHROM <- as.character(gwas_data$CHROM)
gwas_data$CHROM <- sub("^0+(?!$)", "", gwas_data$CHROM, perl = TRUE)

vcf_header <- c(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
)

vcf_body <- gwas_data[, .(CHROM, GENPOS, ID = ".", ALLELE0, ALLELE1, QUAL = ".", FILTER = ".", INFO = paste0("AF=", A1FREQ), FORMAT = ".")]

writeLines(vcf_header, con = paste0(var_path, "/all_traits_in_ld_VEP_input.vcf"))
write.table(vcf_body, paste0(var_path, "/all_traits_in_ld_VEP_input.vcf"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

vcf_body_minimal <- gwas_data[, .(CHROM, GENPOS, ID = ".", ALLELE0, ALLELE1)]
setnames(vcf_body_minimal, c("GENPOS", "ALLELE0", "ALLELE1"), c("POS", "REF", "ALT"))
vcf_body_minimal <- unique(vcf_body_minimal)
write.table(vcf_body_minimal, file = paste0(var_path, "/all_traits_in_ld_CADD_input.vcf"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
