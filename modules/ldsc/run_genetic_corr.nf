/*
 * modules/ldsc/run_genetic_corr.nf
 *
 * Purpose:
 *   Compute pairwise genetic correlations across traits using LDSC outputs
 *   and produce a consolidated log used by higher-level subworkflows.
 *
 * Inputs:
 *   - tuples describing per-trait h2 outputs and filepaths
 *
 * Outputs:
 *   - per-trait genetic correlation log files
 */
nextflow.enable.dsl=2

process RUN_GENETIC_CORR {
  label 'medium'
  tag { trait }

  publishDir "${params.output_path}/Loci_Preprocessing/LDSC_results",
             mode: 'copy',
             overwrite: true,
             pattern: "*genetic_correlation_ldsc_res_ukb*"

  input:
  tuple val(trait), path(runroot), val(trait_file), val(other_trait_files)

  output:
  tuple val(trait), path(runroot), path("LDSC_results/genetic_correlation_ldsc_res_ukb_${trait}.log")

  script:
  """
  set -euo pipefail

  WORK_DIR=\$PWD
  cd "${runroot}"

  mkdir -p "${params.output_path}/Loci_Preprocessing/LDSC_results"
  source config_shell.sh

  # Use explicit Python 2.7 interpreter from config
  ${params.py27_python} "\${LDSC_PATH}" \
    --rg "${trait_file},${other_trait_files}" \
    --ref-ld "\${LDSC_UKBB_DATA}" \
    --w-ld   "\${LDSC_UKBB_DATA}" \
  --out "${params.output_path}/Loci_Preprocessing/LDSC_results/genetic_correlation_ldsc_res_ukb_${trait}"

  mkdir -p "\$WORK_DIR/LDSC_results"
  cp "${params.output_path}/Loci_Preprocessing/LDSC_results/genetic_correlation_ldsc_res_ukb_${trait}"* "\$WORK_DIR/LDSC_results/"
  """
}
