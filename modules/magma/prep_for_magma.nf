/*
 * modules/magma/prep_for_magma.nf
 *
 * Purpose:
 *   Prepare SNP RSID lists for MAGMA by running the R helper that maps
 *   GWAS variants to RSIDs and writes trait-specific *_rsid_magma.txt files.
 *
 * Inputs:
 *   - tuple runroot_path, trait, rsids_file, _is_mtag
 *
 * Outputs:
 *   - ${trait}_*_rsid_magma.txt
 */
nextflow.enable.dsl=2

process PREP_FOR_MAGMA {
  label 'light'
  tag "${trait}"

  publishDir "${params.output_path}/MAGMA",
             mode: 'copy',
             overwrite: true,
             pattern: "*_rsid_magma.txt"

  input:
  tuple val(runroot_path), val(trait), path(rsids_file), val(_is_mtag)

  output:
  tuple val(runroot_path), val(trait), val(false), path("${trait}_*_rsid_magma.txt")

  script:
  def script_path = "${workflow.projectDir}/bin/genes/prep_SNPs_for_magma.R"
  // compute basename in Groovy to avoid embedding shell $(...) which breaks Groovy string parsing
  def rsids_basename = rsids_file.getName()
  """
  set -euo pipefail

  WORK_DIR=\$PWD
  # Copy the staged rsids file into the runroot so R scripts can access it
  cp "${rsids_file}" "${runroot_path}/" || true

  cd ${runroot_path}

  mkdir -p ${params.output_path}/MAGMA
  
  # Use the shared config_R.R that was created by CREATE_CONFIG_SHIMS
  # It already contains gwas_input_path and magma_output_path
  # Pass the rsids filename as a third argument so the R script can locate the input
  Rscript ${script_path} "${trait}" "false" "${rsids_basename}"

  # Stage outputs back to work dir
  cp ${params.output_path}/MAGMA/${trait}_*_rsid_magma.txt "\$WORK_DIR/" || true
  """
}
