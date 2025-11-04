/*
 * modules/prep/initial_qc.nf
 *
 * Purpose:
 *   Run initial QC scripts on the merged per-trait GWAS files and publish
 *   the prepped file used by downstream liftover/rsid steps.
 *
 * Inputs:
 *   - tuple(trait, runroot, script_initial_qc)
 *
 * Outputs:
 *   - <trait>_regenie_allchr.txt (published to GWAS_Preprocessing)
 */
nextflow.enable.dsl=2
process INITIAL_QC {
  label 'light'
  label 'prep_batch5'
  tag { trait }

  input:
    tuple val(trait), path(runroot), path(script_initial_qc)

  output:
  tuple val(trait), path(runroot), path("runroot/${trait}_regenie_allchr.txt")

  script:
  """
  set -euo pipefail
  SCRIPT_PATH="\$(readlink -f "$script_initial_qc" 2>/dev/null || realpath "$script_initial_qc")"
  cd ${runroot}
      "\$SCRIPT_PATH" \
         --in "${trait}_regenie_allchr.txt" \
         --out "${params.output_path}/GWAS_Preprocessing/${trait}_regenie_allchr.txt"
  ln -sf "${params.output_path}/GWAS_Preprocessing/${trait}_regenie_allchr.txt" "./${trait}_regenie_allchr.txt"
  """

  stub:
  """
  mkdir -p runroot
  : > runroot/${trait}_regenie_allchr.txt
  """
}
