/*
 * modules/loci_groups/identify_close_pvals.nf
 *
 * Purpose:
 *   Identify lead SNPs that have close p-values in other traits and produce
 *   a CSV of loci blocks with close p-values for downstream aggregation.
 *
 * Inputs:
 *   - tuple(runroot_path, loci_blocks)
 *
 * Outputs:
 *   - All_loci_blocks_with_close_pvals.csv (published to Loci_Preprocessing)
 */
nextflow.enable.dsl=2

process IDENTIFY_CLOSE_PVALS {
  label 'medium'
  tag "all_traits"

  publishDir "${params.output_path}/Loci_Preprocessing",
             mode: 'copy',
             overwrite: true,
             pattern: "All_loci_blocks_with_close_pvals.csv"

  input:
  tuple val(runroot_path), path(loci_blocks)

  output:
  path "All_loci_blocks_with_close_pvals.csv"

  script:
  def script_path = "${workflow.projectDir}/bin/loci/identify_leads_with_close_pvals.R"
  
  """
  set -euo pipefail
  
  # Store the work directory
  WORK_DIR=\$PWD
  
  # Run the script from within the runroot
  cd ${runroot_path}
  
  # Use the shared config_R.R that was created by CREATE_CONFIG_SHIMS
  # It already contains all required variables including all_loci, traits, gwas_path
  
  # Run the close p-values identification script
  Rscript ${script_path}
  
  # The script modifies all_loci in place, so copy the result
  cp "${params.output_path}/Loci_Preprocessing/Singletrait_all_loci.csv" "\$WORK_DIR/All_loci_blocks_with_close_pvals.csv"
  """
}
