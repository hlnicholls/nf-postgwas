/*
 * modules/pops/clean_results.nf
 *
 * Purpose:
 *   Clean and aggregate PoPS .preds files into a unified features-cleaned
 *   CSV used by downstream per-locus aggregation steps.
 *
 * Inputs:
 *   - tuple(runroot_path, pops_preds)
 *
 * Outputs:
 *   - *pops_results_all_features_cleaned.csv (published to Pops)
 */
nextflow.enable.dsl=2

process CLEAN_POPS_RESULTS {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Pops",
             mode: 'copy',
             overwrite: true,
             pattern: "*pops_results_all_features_cleaned.csv"

  input:
  tuple val(runroot_path), path(pops_preds)

  output:
  tuple val(runroot_path), path("*pops_results_all_features_cleaned.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/genes/clean_pops_results.R"
  """
  set -euo pipefail

  WORK_DIR=\$PWD
  cd ${runroot_path}

  # Use shared config (created by CREATE_CONFIG_SHIMS)
  Rscript ${script_path}

  # Stage cleaned outputs
  cp ${params.output_path}/Pops/*pops_results_all_features_cleaned.csv "\$WORK_DIR/" || true
  """
}
