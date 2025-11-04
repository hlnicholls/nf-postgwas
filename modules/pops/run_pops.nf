/*
 * modules/pops/run_pops.nf
 *
 * Purpose:
 *   Run PoPS per trait to generate gene-level priority scores using MAGMA
 *   results and PoPS features. Publishes *.preds files to the Pops output.
 *
 * Inputs:
 *   - tuple(runroot_path, trait, _is_mtag)
 *
 * Outputs:
 *   - <trait>_pops_all_features.preds (published to Pops)
 */
nextflow.enable.dsl=2

process RUN_POPS {
  label 'heavy'
  tag "${trait}"

  publishDir "${params.output_path}/Pops",
             mode: 'copy',
             overwrite: true,
             pattern: "*.preds"

  input:
  tuple val(runroot_path), val(trait), val(_is_mtag)

  output:
  tuple val(runroot_path), val(trait), path("*.preds")

  script:
  // Fixed, non-configurable base paths
  def outBase       = (params.output_path?.startsWith('/')) ? params.output_path
                                                           : "${runroot_path}/${params.output_path ?: 'results'}"
  def pops_root     = "${params.databases}/pops"
  def pops_py       = "python3"
  def pops_script   = "${pops_root}/pops.py"
  def pops_annot    = "${pops_root}/example/data/utils/gene_annot_jun10.txt"
  def pops_features = "${pops_root}/features_munged_all/pops_features_all"
  def pops_control  = "${pops_root}/example/data/utils/features_jul17_control.txt"
  def n_chunks      = 2

  def popsOutDir  = "${outBase}/Pops"
  def magmaResDir = "${outBase}/MAGMA/Results"
  def magmaPrefix = "magma_${trait}"
  def outPrefix   = "${trait}_pops_all_features"

  """
  set -euo pipefail

  WORK_DIR="\$PWD"
  cd "${pops_root}"

  mkdir -p "${popsOutDir}"

  echo "Processing phenotype: ${trait}"
  echo "Using PoPS at: \$(pwd)"
  echo "MAGMA results prefix: ${magmaResDir}/${magmaPrefix}"
  echo "Output prefix: ${popsOutDir}/${outPrefix}"

  "${pops_py}" "${pops_script}" \\
    --gene_annot_path "${pops_annot}" \\
    --feature_mat_prefix "${pops_features}" \\
    --num_feature_chunks "${n_chunks}" \\
    --magma_prefix "${magmaResDir}/${magmaPrefix}" \\
    --control_features_path "${pops_control}" \\
    --out_prefix "${popsOutDir}/${outPrefix}"

  cp "${popsOutDir}/${outPrefix}.preds" "\${WORK_DIR}/"
  """
}
