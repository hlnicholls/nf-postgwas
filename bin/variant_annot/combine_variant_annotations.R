library(dplyr)
library(data.table)


source("config_R.R")

# Find all files in the input directory
all_files <- list.files(var_path, full.names = TRUE)

# Exclude "all_traits_in_ld_VEP_revel_annotated.txt" and filter valid file types
files_to_merge <- all_files[!grepl("all_traits_in_ld_VEP_revel_annotated.txt|all_variant_annotations.csv", all_files) &
                              (grepl("\\.csv$", all_files) | grepl("\\.tsv$", all_files) | grepl("\\.txt$", all_files) | grepl("\\.gz$", all_files))]
# Initialize an empty DataFrame for merging
merged_df <- NULL

# Read and merge all files based on the 'SNP' column
for (file_path in files_to_merge) {
  # Determine file type and read accordingly
  if (grepl("\\.gz$", file_path)) {
    df <- fread(cmd = paste("zcat", file_path))
  } else if (grepl("\\.csv$", file_path)) {
    df <- fread(file_path)
  } else if (grepl("\\.tsv$", file_path) | grepl("\\.txt$", file_path)) {
    df <- fread(file_path, sep = "\t")
  } else {
    stop(paste("Unsupported file format:", file_path))
  }
  print('Merging:')
  print(file_path)
  #print(head(df))
  # Merge files on 'SNP' column, keeping all SNPs
  if (is.null(merged_df)) {
    merged_df <- df
  } else {
    merged_df <- merge(merged_df, df, by = "SNP", all=TRUE)  # Use full_join to keep all rows
  }
}

merged_df <- merged_df %>%
  distinct()
  
#print('Merged Annotations:')
#print(head(merged_df, 10))

# Read the loci file to be annotated
loci_df <- fread(var_file, sep = "\t")

# Merge loci DataFrame with the merged annotation DataFrame, keeping only SNPs from loci_df
final_merged_df <- merge(loci_df, merged_df, by = "SNP", all.x=TRUE)

# Collapse and clean duplicate values
final_merged_df <- final_merged_df %>%
  group_by(SNP, Phenotype, lead_snp, Method) %>%
  summarise_all(~paste(unique(na.omit(.)), collapse = ","))

# Print the final merged DataFrame after collapsing and cleaning
print('Final Merged DataFrame with collapsed and cleaned values:')
print(head(final_merged_df, 10))
print(paste('Shape of original loci dataframe:', dim(loci_df)[1], dim(loci_df)[2]))
print(paste('Shape of final merged data - after collapsing duplicate loci (if from multiple phenotypes and single and multi-trait):', dim(final_merged_df)[1], dim(final_merged_df)[2]))
# Order final_merged_df by CHROM and lead_snp
final_merged_df <- final_merged_df[order(final_merged_df$CHROM, final_merged_df$lead_snp), ]

data.table::fwrite(final_merged_df, file.path(var_path, "all_variant_annotations.csv"))