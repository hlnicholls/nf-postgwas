# ============================================================
# GLOBAL SETTINGS
# ============================================================
traits <- c('smoketrait')
sample_sizes <- c(123)
home_path <- here::here()
log_path <- file.path(here::here(), "Logs")
output_path <- "/postgwas/results/postGWAS_Output"
databases <- "/postgwas/databases/"

# ============================================================
# SECTION 1: INPUT DIRECTORIES
# ============================================================
# Processed GWAS summary stats from PREP_FLOW
gwas_input_path <- file.path(output_path, "GWAS_Preprocessing")

# ============================================================
# SECTION 2: OUTPUT DIRECTORIES
# ============================================================
loci_output_path         <- file.path(output_path, "Loci_Preprocessing")
plots_output_path        <- file.path(output_path, "GWAS_Plots")
coloc_output_path        <- file.path(output_path, "Colocalisation")
enrichment_output_path   <- file.path(output_path, "Enrichment")
magma_output_path        <- file.path(output_path, "MAGMA")
pops_output_path         <- file.path(output_path, "Pops")
depict_output_path       <- file.path(output_path, "DEPICT")
prs_output_path          <- file.path(output_path, "PRS")
phewas_output_path       <- file.path(output_path, "PheWAS")
prioritised_genes_path   <- file.path(output_path, "Prioritisation")
druggability_output_path <- file.path(output_path, "Druggability")
finemapping_path         <- file.path(output_path, "Fine-mapping/Wakefield")
var_path                 <- file.path(output_path, "Variant_annotation")
annotated_credible_set_path <- file.path(output_path, "Credible_set_annotation")
pleiotropy_outpath       <- file.path(output_path, "Pleiotropy")

# ============================================================
# SECTION 3: SPECIFIC FILE PATHS
# ============================================================
all_loci    <- file.path(loci_output_path, "Singletrait_all_loci.csv")
# Canonical reference used downstream
var_file    <- file.path(loci_output_path, "all_traits_loci_38_with_ld_genes.txt")

locus_blocks_output <- file.path(loci_output_path, "All_loci_blocks.csv")
ldsc_results_path   <- file.path(loci_output_path, "LDSC_results")

# Prioritisation / coloc / FM
prioritised_genes  <- file.path(prioritised_genes_path, "Prioritised_genes.csv")
coloc_all_variants <- file.path(coloc_output_path, "HF_coloc_all_variants_pp4.csv")
fm_vars            <- file.path(loci_output_path, "all_traits_FM_vars_annotated.txt")
annotated_vars     <- file.path(loci_output_path, "all_traits_vars_annot_38.txt")
HiC_path           <- file.path(databases, "hic/All_HiC_enhancer_gene.txt")

# Enrichment
gprofiler_results <- file.path(enrichment_output_path, "gprofiler_results_table.csv")
gene_targets      <- file.path(databases, "opentargets/Gene_target_prioritisation.csv")
impc_results_path <- file.path(enrichment_output_path, "IMPC_results_table.csv")

# Fine-mapping sets
wakefield_cs <- file.path(finemapping_path, "All_traits_credible_sets.txt")
susie_cs     <- file.path(output_path, "Fine-mapping/SuSiE/All_traits_credible_sets.txt")
carma_cs     <- file.path(output_path, "Fine-mapping/CARMA/All_traits_credible_sets.txt")
credidble_sets_path <- file.path(annotated_credible_set_path, "Annotated_CredibleSets_and_LD08.csv")

# DEPICT
depict_path         <- depict_output_path
depict_results_path <- depict_output_path
depict_cfg_template <- file.path(databases, "depict/config_phenotype_1e5.cfg")

# PRS
PRS_path       <- file.path(prs_output_path, "PRScs")
PRS_plink_PATH <- file.path(prs_output_path, "PRS_plink_scores")

# PheWAS
prs_cs_score_path <- PRS_plink_PATH
phewas_results    <- phewas_output_path

# ============================================================
# BACKWARD COMPATIBILITY ALIASES
# ============================================================
gwas_path         <- gwas_input_path
plot_outpath      <- plots_output_path
coloc_path        <- coloc_output_path
enrichment_path   <- enrichment_output_path
magma_path        <- magma_output_path
pops_results_path <- pops_output_path

# ============================================================
# SPECIALIZED PATHS
# ============================================================
# Plotting
locuszoom_path <- file.path(plots_output_path, "Regional_plots")

# Colocalisation paths & GTEx
HF_gwas           <- file.path(databases, "coloc_gwas/Levin_HF_lifted.txt")
DCM_gwas          <- file.path(databases, "coloc_gwas/DCM_gwas_lifted.txt")
custom_gwas       <- ""
custom_trait_name <- "MyTrait"
eqtl_in           <- file.path(coloc_output_path, "eQTL")
gtex_path         <- file.path(databases, "gtex_v8_data/")
eqtl_dir          <- file.path(databases, "gtex_v8_data/")
tissues <- c("Artery_Aorta", "Artery_Coronary", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
HiC_tissues <- c()
genetic_map_hg38 <- "/postgwas/databases/genetic_map/genetic_map_hg38_withX_broad_eagle.txt"
disease_terms <- c()
GTEx_tissue_terms  <- c()
disease_regex <- paste(disease_terms, collapse="|")
tissue_regex  <- paste(GTEx_tissue_terms,  collapse="|")
