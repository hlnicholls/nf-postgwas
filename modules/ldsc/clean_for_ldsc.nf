/*
 * modules/ldsc/clean_for_ldsc.nf
 *
 * Purpose:
 *   Prepare and clean GWAS summary statistics so they are compatible with
 *   LDSC input format. This step runs per trait and emits cleaned files.
 *
 * Inputs:
 *   - per-trait triples (trait, runroot, rsids_txt)
 *
 * Outputs:
 *   - cleaned files for LDSC heritability and genetic correlation
 */
nextflow.enable.dsl=2


process CLEAN_FOR_LDSC {
  label 'light'
  tag { trait }

  publishDir "${params.output_path}/GWAS_Preprocessing",
             mode: 'copy',
             overwrite: true,
             pattern: "*GWAS_37_corr*.txt"

  input:
  tuple val(trait), path(runroot), path(rsids_file), path(script_clean_for_ldsc)

  output:
  tuple val(trait), path(runroot), path("${trait}_GWAS_37_corr.txt")

  script:
  """
  set -euo pipefail

  # Store work directory
  WORK_DIR=\$PWD

  # Resolve script path before changing directories
  SCRIPT_PATH="\$(readlink -f "$script_clean_for_ldsc" 2>/dev/null || realpath "$script_clean_for_ldsc")"

  # Run the LDSC cleaning script from within the runroot
  cd ${runroot}
    "\$SCRIPT_PATH"

  # Copy output file back to work directory
  cp ${params.output_path}/GWAS_Preprocessing/${trait}_GWAS_37_corr.txt "\$WORK_DIR/"
  """
}
