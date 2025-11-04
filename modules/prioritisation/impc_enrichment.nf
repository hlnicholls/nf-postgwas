/*
 * modules/prioritisation/impc_enrichment.nf
 *
 * Purpose:
 *   Run IMPC enrichment analyses on the prioritised gene list and publish
 *   ranked tables and summary plots.
 *
 * Inputs:
 *   - tuple runroot_path, prioritised_genes_csv
 *
 * Outputs:
 *   - IMPC_genes_ranked_enriched.csv
 *   - IMPC_results_table.csv
 *   - Genes_ranked_IMPC_enrichment.png
 */
nextflow.enable.dsl=2

process IMPC_ENRICHMENT {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Enrichment",
             mode: 'copy',
             overwrite: true

  input:
  // From GENE_PRIORITISATION: (runroot_path, Prioritised_genes.csv)
  tuple val(runroot_path), path(prioritised_genes_csv)

  output:
  // Emit all three files as a single tuple channel
  path "IMPC_genes_ranked_enriched.csv"
  path "IMPC_results_table.csv"
  path "Genes_ranked_IMPC_enrichment.png"

  script:
  def script_path = "${workflow.projectDir}/bin/genes/IMPC_enrichment.R"
  """
  set -euo pipefail

  WORK_DIR="\$PWD"

  # Validate inputs
  [ -s "${prioritised_genes_csv}" ] || { echo "Missing prioritised genes CSV: ${prioritised_genes_csv}" >&2; exit 1; }
  [ -d "${runroot_path}" ]          || { echo "Missing runroot directory: ${runroot_path}" >&2; exit 1; }

  cd "${runroot_path}"

  # Ensure enrichment directory exists (script writes there using config)
  mkdir -p "${params.output_path}/Enrichment"

  echo "IMPC_ENRICHMENT: running in \$(pwd)"
  echo "IMPC_ENRICHMENT: invoking R script: ${script_path}"

  Rscript "${script_path}"

  # Stage outputs back to work dir so Nextflow can capture them
  cp "${params.output_path}/Enrichment/IMPC_genes_ranked_enriched.csv"       "\$WORK_DIR/" || true
  cp "${params.output_path}/Enrichment/IMPC_results_table.csv"               "\$WORK_DIR/" || true
  cp "${params.output_path}/Enrichment/Genes_ranked_IMPC_enrichment.png"     "\$WORK_DIR/" || true

  # Final sanity check
  ls -l "\$WORK_DIR"/IMPC_* "\$WORK_DIR"/Genes_ranked_IMPC_enrichment.png
  """
}
