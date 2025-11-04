/*
 * modules/variant_annot/combine_annotations.nf
 *
 * Purpose:
 *   Combine per-annotation outputs (CADD, pCHiC, RegulomeDB, VEP) into a
 *   single consolidated CSV for downstream use and publish it to
 *   Variant_annotation/all_variant_annotations.csv.
 *
 * Inputs:
 *   - val runroot_path
 *
 * Outputs:
 *   - all_variant_annotations.csv
 */
nextflow.enable.dsl=2

process VARIANT_ANNOT_COMBINE {
  label 'medium'
  tag "all_traits"

  publishDir "${params.output_path}/Variant_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "all_variant_annotations.csv"

  input:
  val runroot_path

  output:
  path "all_variant_annotations.csv"

  script:
  def script_path = "${workflow.projectDir}/bin/variant_annot/combine_variant_annotations.R"
  def runroot = runroot_path  // Store in a variable for interpolation
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  cd ${runroot}

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Variant_annotation"

  # Use shared config (created by CREATE_CONFIG_SHIMS)
  Rscript ${script_path}

  # Stage output back to work dir
  cp "${params.output_path}/Variant_annotation/all_variant_annotations.csv" "\$WORK_DIR/"
  """
}
