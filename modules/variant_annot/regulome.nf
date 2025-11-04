/*
 * modules/variant_annot/regulome.nf
 *
 * Purpose:
 *   Annotate variants with RegulomeDB scores via an R helper and publish
 *   the combined results to Variant_annotation.
 *
 * Inputs:
 *   - tuple(runroot_path, ld_genes_file)
 *
 * Outputs:
 *   - all_traits_in_ld_regulomedb.txt
 */
nextflow.enable.dsl=2

process VARIANT_ANNOT_REGULOME {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Variant_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "all_traits_in_ld_regulomedb.txt"

  input:
  tuple val(runroot_path), path(ld_genes_file)

  output:
  tuple val(runroot_path), path("all_traits_in_ld_regulomedb.txt")

  script:
  def script_path = "${workflow.projectDir}/bin/variant_annot/regulomeDB_score.R"
  def outputFile = file("${params.output_path}/Variant_annotation/all_traits_in_ld_regulomedb.txt")
  
  if (outputFile.exists() && outputFile.size() > 0) {
    """
    set -euo pipefail
    
    echo "⊘ RegulomeDB results already exist, copying existing file..."
    cp "${outputFile}" all_traits_in_ld_regulomedb.txt
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
  cp "${params.output_path}/Variant_annotation/all_traits_in_ld_regulomedb.txt" "\$WORK_DIR/"
    """
  }
}
