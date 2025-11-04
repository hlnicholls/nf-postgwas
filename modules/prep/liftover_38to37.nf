/*
 * modules/prep/liftover_38to37.nf
 *
 * Purpose:
 *   Run the liftover helper to convert GWAS coordinates from GRCh38 to GRCh37
 *   for downstream compatibility with older reference resources.
 *
 * Inputs:
 *   - tuple(trait, runroot, qc_merged, script_liftover)
 *
 * Outputs:
 *   - <trait>_38_37.txt (symlinked into runroot)
 */
nextflow.enable.dsl=2
process LIFTOVER_38_TO_37 {
  label 'light'
  label 'prep_batch5'
  tag { trait }

  input:
  tuple val(trait), path(runroot), path(qc_merged), path(script_liftover)

  output:
  tuple val(trait), path(runroot), path("runroot/${trait}_38_37.txt")

  script:
  """
  set -euo pipefail
  SCRIPT_PATH="\$(readlink -f "$script_liftover" 2>/dev/null || realpath "$script_liftover")"
  cd ${runroot}
    "\$SCRIPT_PATH" \
      --in "${qc_merged}" \
      --out "${params.output_path}/GWAS_Preprocessing/${trait}_38_37.txt"
  ln -sf "${params.output_path}/GWAS_Preprocessing/${trait}_38_37.txt" "./${trait}_38_37.txt"
  """

  stub:
  """
  mkdir -p runroot
  : > runroot/${trait}_38_37.txt
  """
}
