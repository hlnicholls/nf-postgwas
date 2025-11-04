/*
 * modules/utils/make_config_shims.nf
 *
 * Purpose:
 *   Create per-trait runroot configuration shims for R/Python/shell helper
 *   scripts. Writes config_R.R, config_python.py and config_shell.sh into
 *   the runroot directory so downstream scripts can locate paths.
 *
 * Inputs:
 *   - tuple trait, sample_size, genome_build
 *
 * Outputs:
 *   - runroot (directory) containing config files
 */
nextflow.enable.dsl=2

process CREATE_CONFIG_SHIMS {
  label 'light'
  tag { trait }

  input:
  tuple val(trait), val(sample_size), val(genome_build)

  output:
  tuple val(trait), path('runroot', type: 'dir')

  script:
  // --- Stable prefixes & values resolved by Groovy (Nextflow) at compile time ---
  def dbPrefix        = params.databases?.endsWith('/') ? params.databases : "${params.databases}/"
  def geneticMapPath  = params.genetic_map_hg38 ?: "${dbPrefix}genetic_map/genetic_map_hg38_withX_broad_eagle.txt"
  def gtexTissuesR    = (params.gtex_tissues ?: []).collect { "\"${it}\"" }.join(', ')
  def hicTissuesR     = (params.HiC_tissues   ?: []).collect { "\"${it}\"" }.join(', ')

  // Render user-provided disease/tissue terms for R and Python shims
  def diseaseTermsR   = (params.disease_terms ?: []).collect { "\"${it}\"" }.join(', ')
  def tissueTermsR    = (params.GTEx_tissue_terms  ?: []).collect { "\"${it}\"" }.join(', ')
  def pyDisease       = (params.disease_terms ?: []).collect { "\"${it}\"" }.join(', ')
  def pyTissue        = (params.GTEx_tissue_terms  ?: []).collect { "\"${it}\"" }.join(', ')

  """
  set -euo pipefail

  # Create runroot structure
  mkdir -p runroot/Logs/GWAS_Preprocessing
  mkdir -p "${params.output_path}/GWAS_Preprocessing"

  # -------------------------------
  # Write R config (literal script)
  # -------------------------------
  cat > runroot/config_R.R << 'EOF'
# ============================================================
# GLOBAL SETTINGS
# ============================================================
traits <- c('${trait}')
sample_sizes <- c(${sample_size})
home_path <- here::here()
log_path <- file.path(here::here(), "Logs")
output_path <- "${params.output_path}"
databases <- "${dbPrefix}"

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
custom_gwas       <- "${params.custom_gwas_path}"
custom_trait_name <- "${params.custom_trait_name}"
eqtl_in           <- file.path(coloc_output_path, "eQTL")
gtex_path         <- file.path(databases, "gtex_v8_data/")
eqtl_dir          <- file.path(databases, "gtex_v8_data/")
EOF

  # Append dynamic R bits that depend on params lists/paths
  echo 'tissues <- c(${gtexTissuesR})'            >> runroot/config_R.R
  echo 'HiC_tissues <- c(${hicTissuesR})'         >> runroot/config_R.R
  echo 'genetic_map_hg38 <- "${geneticMapPath}"'  >> runroot/config_R.R
  # User-provided keyword sets & compiled regex (for grepl in R)
  echo 'disease_terms <- c(${diseaseTermsR})'     >> runroot/config_R.R
  echo 'GTEx_tissue_terms  <- c(${tissueTermsR})'      >> runroot/config_R.R
  echo 'disease_regex <- paste(disease_terms, collapse="|")' >> runroot/config_R.R
  echo 'tissue_regex  <- paste(GTEx_tissue_terms,  collapse="|")' >> runroot/config_R.R

  # -------------------------------
  # Write Python config (literal) — VALID PYTHON ONLY
  # -------------------------------
  cat > runroot/config_python.py << EOF
traits = ["${trait}"]
sample_sizes = [${sample_size}]
genome_build = "${genome_build}"

home_path = "."
gwas_input_dir = "${params.gwas_input_dir}"
# legacy Step1_initial_QC.py compatibility
input_path = "${params.output_path}/GWAS_Preprocessing"
input_suffix = "_assoc_regenie_allchr.txt"

output_path = "${params.output_path}"
databases = "${dbPrefix}"
log_path = f"{home_path}/Logs"

# LD and loci paths for collate/annotate scripts
loci_preprocessing = f"{output_path}/Loci_Preprocessing"
ld_path = f"{loci_preprocessing}/Single_trait_LD"
loci_path = f"{loci_preprocessing}/Singletrait_all_loci.csv"

# Variant annotation and fine-mapping paths
var_path = f"{output_path}/Variant_annotation"
finemapping_path = f"{output_path}/Fine-mapping/Wakefield"
annotated_credible_set_path = f"{output_path}/Credible_set_annotation"
pleiotropy_outpath = f"{output_path}/Pleiotropy/"

# PRS and PheWAS paths
PRS_path = f"{output_path}/PRS/PRScs"
PRS_plink_PATH = f"{output_path}/PRS/PRS_plink_scores"
mtag_results_qc_path = None  # set to MTAG path if using MTAG results

# Fine-mapping specific files
wakefield_cs = f"{finemapping_path}/All_traits_credible_sets.txt"
susie_cs = f"{output_path}/Fine-mapping/SuSiE/All_traits_credible_sets.txt"
carma_cs = f"{output_path}/Fine-mapping/CARMA/All_traits_credible_sets.txt"

# Optional keyword sets also available in Python (not used by collate_LD.py)
disease_terms = [${pyDisease}]
GTEx_tissue_terms  = [${pyTissue}]
disease_regex = "|".join(disease_terms) if disease_terms else ""
tissue_regex  = "|".join(GTEx_tissue_terms)  if GTEx_tissue_terms  else ""
EOF

  # ---------------------------------------------------
  # Write shell config with REAL values (interpolated)
  # Use SINGLE-QUOTED heredoc to prevent premature expansion.
  # ---------------------------------------------------
  cat > runroot/config_shell.sh << 'EOF'
# =============================
# Shell configuration shim
# =============================

# Traits and sizes
TRAITS=("${trait}")
SAMPLE_SIZES="${sample_size}"        # vector-friendly string for your scripts
SAMPLE_SIZE="${sample_size}"         # some scripts expect a singular

# Core paths
HOME_PATH="."
OUTPUT_PATH="${params.output_path}"
COMBINED_REGENIE_PATH="\${OUTPUT_PATH}/GWAS_Preprocessing"
REGENIE_FOLDER="${params.gwas_input_dir}"
DATABASES="${dbPrefix}"              # with trailing slash
DATABASES_DIR="${dbPrefix}"          # alias used by legacy scripts
LOG_PATH="\${HOME_PATH}/Logs"
ERROR_PATH="\${HOME_PATH}/Errors"

# Optional: reference panel / downsample toggle
DOWNSAMPLE="${params.ld_reference_panel}"
export DOWNSAMPLE


# Software
PLINK1_PATH="plink"
PLINK2_PATH="plink2"
LDSC_PATH="\${DATABASES_DIR}ldsc/ldsc.py"
export PLINK1_PATH PLINK2_PATH

# -----------------------------
# Step 1 - Loci Preprocessing
# -----------------------------
LOCI_PREPROCESSING="\${OUTPUT_PATH}/Loci_Preprocessing"

# Single-trait outputs
SINGLE_TRAIT_LOCI="\${LOCI_PREPROCESSING}/Singletrait_all_loci.csv"
SINGLE_TRAIT_LOCI_HCM="\${LOCI_PREPROCESSING}/Single_trait_LD/missing_loci.csv"
SINGLE_TRAIT_LD="\${LOCI_PREPROCESSING}/Single_trait_LD"

# LDSC folders & data
LDSC_OUT="\${LOCI_PREPROCESSING}/LDSC_results"
LDSC_UKBB_DATA="\${DATABASES_DIR}ldsc/UKBB.ALL.ldscore/UKBB.EUR"

# LD output per locus
ALL_TRAIT_LD="\${LOCI_PREPROCESSING}/All_trait_LD"

# Lead SNPs lists (single vs MTAG)
ALL_TRAIT_LOCI="\${LOCI_PREPROCESSING}/Single_and_MTAG_all_loci.csv"
ALL_TRAIT_UNIQUE_LOCI="\${LOCI_PREPROCESSING}/Single_trait_loci_unique_across_all_traits.csv"

# Variant annotation integration
VEP_OUTPUT_NAME="all_traits_in_ld_VEP_revel_annotated.txt"
VAR_PATH="\${OUTPUT_PATH}/Variant_annotation"

# Fine-mapping
LD_MATRICES="\${OUTPUT_PATH}/Fine-mapping/LD_matrices"

# MAGMA
MAGMA_1000G="\${DATABASES_DIR}MAGMA/1000G.EUR"
MAGMA_annot="\${DATABASES_DIR}MAGMA/magma_0kb.genes.annot"
MAGMA_PATH="\${OUTPUT_PATH}/MAGMA"

# PoPS
POPS_SCRIPT="\${DATABASES_DIR}pops/pops.py"
POPS_ANNOT="\${DATABASES_DIR}pops/example/data/utils/gene_annot_jun10.txt"
POPS_FEATURES_ALL="\${DATABASES_DIR}pops/features_munged_all/pops_features_all"
POPS_CONTROL="\${DATABASES_DIR}pops/example/data/utils/features_jul17_control.txt"
POPS_PATH="\${OUTPUT_PATH}/Pops"

# DEPICT
DEPICT_PATH="\${OUTPUT_PATH}/DEPICT"

# PRS (two layouts preserved for compatibility)
PRS_SCRIPT="\${DATABASES_DIR}PRScs-master/PRScs.py"
PRS_PATH="\${OUTPUT_PATH}/PRS/PRScs"
PRS_plink_PATH="\${OUTPUT_PATH}/PRS/PRS_plink_scores"

# Newer PRS paths used by R config too
PRS_PATH_NEW="\${OUTPUT_PATH}/PRS/PRScs"
PRS_PLINK_PATH_NEW="\${OUTPUT_PATH}/PRS/PRS_plink_scores"
EOF
  """

  stub:
  """
  mkdir -p runroot
  echo "Stub: CREATE_CONFIG_SHIMS for ${trait}"
  """
}
