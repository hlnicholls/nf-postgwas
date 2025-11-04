/*
 * modules/variant_annot/pchic.nf
 *
 * Purpose:
 *   Run promoter-capture Hi-C (pCHiC) annotation helper and publish the
 *   aggregated per-trait HiC annotation TSV to Variant_annotation.
 *
 * Inputs:
 *   - tuple(runroot_path, ld_genes_file)
 *
 * Outputs:
 *   - all_traits_in_ld_HiC.tsv
 */
nextflow.enable.dsl=2

process VARIANT_ANNOT_PCHIC {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Variant_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "all_traits_in_ld_HiC.tsv"

  input:
  tuple val(runroot_path), path(ld_genes_file)

  output:
  tuple val(runroot_path), path("all_traits_in_ld_HiC.tsv")

  script:
  def script_path = "${workflow.projectDir}/bin/variant_annot/pCHiC_annotation.R"
  def outputFile = file("${params.output_path}/Variant_annotation/all_traits_in_ld_HiC.tsv")
  
  if (outputFile.exists() && outputFile.size() > 0) {
    """
    # pCHiC annotation already exists - copy to work directory
  cp "${params.output_path}/Variant_annotation/all_traits_in_ld_HiC.tsv" .
    """
  } else {
    """
    set -euo pipefail

    WORK_DIR=\$PWD

    # Touch input to ensure it exists
    test -s "${ld_genes_file}" > /dev/null

    cd ${runroot_path}

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Variant_annotation"

    # Use shared config (created by CREATE_CONFIG_SHIMS)
    Rscript ${script_path}

    # Stage output back to work dir
    cp "${params.output_path}/Variant_annotation/all_traits_in_ld_HiC.tsv" "\$WORK_DIR/"
    """
  }
}
