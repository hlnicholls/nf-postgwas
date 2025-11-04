/*
 * modules/variant_annot/cadd_from_source.nf
 *
 * Purpose:
 *   Compute CADD annotations from source R scripts and publish the
 *   combined per-trait CADD CSV to Variant_annotation.
 *
 * Inputs:
 *   - tuple(runroot_path, ld_genes_file)
 *
 * Outputs:
 *   - all_traits_in_ld_CADD.csv
 */
nextflow.enable.dsl=2

process VARIANT_ANNOT_CADD {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Variant_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "all_traits_in_ld_CADD.csv"

  input:
  tuple val(runroot_path), path(ld_genes_file)

  output:
  tuple val(runroot_path), path("all_traits_in_ld_CADD.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/variant_annot/CADD_score_from_source.R"
  def outputFile = "${params.output_path}/Variant_annotation/all_traits_in_ld_CADD.csv"
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Variant_annotation"

  # Check if output already exists to enable skipping
  if [ -f "${outputFile}" ] && [ -s "${outputFile}" ]; then
    echo "⊘ CADD annotation already exists, copying to work directory"
    cp "${outputFile}" "\$WORK_DIR/"
  else
    echo "✓ Running CADD annotation from source"
    
    # Touch input to ensure it exists
    test -s "${ld_genes_file}" > /dev/null
    
    cd ${runroot_path}

    # Use shared config (created by CREATE_CONFIG_SHIMS)
    Rscript ${script_path}

  # Stage output back to work dir
  cp "${outputFile}" "\$WORK_DIR/"
  fi
  """
}
