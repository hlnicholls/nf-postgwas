/*
 * modules/pops/pops_per_loci.nf
 *
 * Purpose:
 *   Aggregate cleaned PoPS predictions per locus using loci block definitions
 *   and LD-annotated gene lists. Publishes top genes and per-locus summaries.
 *
 * Inputs:
 *   - tuple(runroot_path, cleaned_csv)
 *   - tuple(_runroot2, locus_blocks_csv)
 *   - tuple(_runroot3, ld_genes_file)
 *
 * Outputs:
 *   - pops_top_genes_per_locus.txt
 *   - gwas_all_loci_top_pops_genes.txt
 */
nextflow.enable.dsl=2

process POPS_PER_LOCI {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Pops",
             mode: 'copy',
             overwrite: true,
             pattern: "{pops_top_genes_per_locus.txt,gwas_all_loci_top_pops_genes.txt}"

  input:
  // Ensure ordering with upstream steps (from different channels)
  tuple val(runroot_path), path(cleaned_csv)
  tuple val(_runroot2),    path(locus_blocks_csv)
  tuple val(_runroot3),    path(ld_genes_file)

  output:
  tuple path("pops_top_genes_per_locus.txt"), path("gwas_all_loci_top_pops_genes.txt")

  script:
  def script_path = "${workflow.projectDir}/bin/genes/pops_per_loci.R"
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  # Touch inputs to satisfy Nextflow 'unused variable' checks and ensure presence
  # Do this BEFORE changing directory
  test -s "${cleaned_csv}"       > /dev/null
  test -s "${locus_blocks_csv}" > /dev/null
  test -s "${ld_genes_file}"    > /dev/null

  cd ${runroot_path}

  # Use shared config (created by CREATE_CONFIG_SHIMS)
  Rscript ${script_path}

  # Stage outputs back to work dir
  cp ${params.output_path}/Pops/pops_top_genes_per_locus.txt "\$WORK_DIR/" || true
  cp ${params.output_path}/Pops/gwas_all_loci_top_pops_genes.txt "\$WORK_DIR/" || true
  """
}
