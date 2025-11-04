/*
 * modules/fine_mapping/wakefield.nf
 *
 * Purpose:
 *   Compute Wakefield approximate Bayes factor credible sets per trait
 *   using the supplied variant annotations and publish CSV/TXT outputs.
 *
 * Inputs:
 *   - tuple(runroot_path, variant_annotations)
 *
 * Outputs:
 *   - credible set CSV/TXT files under Fine-mapping/Wakefield
 */
nextflow.enable.dsl=2

process WAKEFIELD_CREDIBLE_SETS {
  label 'medium'
  tag "all_traits"

  publishDir "${params.output_path}/Fine-mapping/Wakefield",
             mode: 'copy',
             overwrite: true,
             pattern: "*.{csv,txt}"

  input:
  tuple val(runroot_path), path(variant_annotations)

  output:
  tuple val(runroot_path), path("*.{csv,txt}")

  script:
  def script_path = "${workflow.projectDir}/bin/fine_mapping/wakefield_credible_sets.R"
  def outputDir = file("${params.output_path}/Fine-mapping/Wakefield")
  
  // Check if Wakefield credible sets results already exist
  def resultsExist = outputDir.exists() && 
                     outputDir.list().any { it.contains('credible_sets') || it.contains('wakefield') }
  
  if (resultsExist) {
    """
  # Wakefield credible sets already exist - copy to work directory
  cp "${params.output_path}/Fine-mapping/Wakefield"/*.csv . 2>/dev/null || true
  cp "${params.output_path}/Fine-mapping/Wakefield"/*.txt . 2>/dev/null || true
    """
  } else {
    """
    set -euo pipefail

    WORK_DIR=\$PWD

    # Touch input to ensure it exists
    test -s "${variant_annotations}" > /dev/null

    cd ${runroot_path}

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Fine-mapping/Wakefield"

    # Use shared config (created by CREATE_CONFIG_SHIMS)
    Rscript ${script_path}

  # Stage outputs back to work dir
  cp "${params.output_path}/Fine-mapping/Wakefield"/*.csv "\$WORK_DIR/" 2>/dev/null || true
  cp "${params.output_path}/Fine-mapping/Wakefield"/*.txt "\$WORK_DIR/" 2>/dev/null || true
    """
  }
}
