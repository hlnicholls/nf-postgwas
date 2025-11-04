/*
 * modules/ldsc/run_heritability.nf
 *
 * Purpose:
 *   Run per-trait heritability estimation using LDSC (LD Score Regression)
 *   and publish logs/results for downstream aggregation.
 *
 * Inputs:
 *   - tuples produced by CLEAN_FOR_LDSC
 *
 * Outputs:
 *   - per-trait LDSC results (published under params.output_path/Loci_Preprocessing/LDSC_results)
 */
nextflow.enable.dsl=2

process RUN_HERITABILITY {
  label 'medium'
  tag { trait }

  publishDir "${params.output_path}/Loci_Preprocessing/LDSC_results",
             mode: 'copy',
             overwrite: true,
             pattern: "*ldsc_res_ukb*"

  input:
  tuple val(trait), path(runroot), path(original_file)

  output:
  tuple val(trait), path(runroot), path("LDSC_results/${trait}_ldsc_res_ukb.log")

  script:
  """
  set -euo pipefail

  # Store work directory
  WORK_DIR=\$PWD

  cd "${runroot}"

  # Ensure LDSC results directory exists
  mkdir -p "${params.output_path}/Loci_Preprocessing/LDSC_results"

  # Source config for environment variables
  source config_shell.sh

  # Use explicit Python 2.7 interpreter from config
  ${params.py27_python} "\${LDSC_PATH}" \
    --h2 "${params.output_path}/GWAS_Preprocessing/${trait}_GWAS_37_corr.txt" \
    --ref-ld "\${LDSC_UKBB_DATA}" \
    --w-ld   "\${LDSC_UKBB_DATA}" \
  --out    "${params.output_path}/Loci_Preprocessing/LDSC_results/${trait}_ldsc_res_ukb"

  # Stage output to match declared Nextflow output
  mkdir -p "\$WORK_DIR/LDSC_results"
  if [ -f "${params.output_path}/Loci_Preprocessing/LDSC_results/${trait}_ldsc_res_ukb.log" ]; then
    cp "${params.output_path}/Loci_Preprocessing/LDSC_results/${trait}_ldsc_res_ukb.log" "\$WORK_DIR/LDSC_results/"
  fi
  """
}
