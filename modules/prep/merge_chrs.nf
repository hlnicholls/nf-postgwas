/*
 * modules/prep/merge_chrs.nf
 *
 * Purpose:
 *   Merge per-chromosome regenie outputs into a single per-trait file
 *   and publish the canonical merged file for later PREP steps.
 *
 * Inputs:
 *   - tuple(trait, runroot, script_merge_chrs)
 *
 * Outputs:
 *   - <trait>_regenie_allchr.txt (symlinked into runroot)
 */
nextflow.enable.dsl=2

process MERGE_CHRS {
  label 'light'
  label 'prep_batch5'
  tag { trait }

  input:
  tuple val(trait), path(runroot), path(script_merge_chrs)

  output:
  tuple val(trait), path("runroot/${trait}_regenie_allchr.txt")

  script:
  """
  set -euo pipefail
    echo "DEBUG: script_merge_chrs path: $script_merge_chrs" >&2
    echo "DEBUG: workflow.projectDir: ${workflow.projectDir}" >&2
  # Resolve absolute path to the script before changing directories
  SCRIPT_PATH="\$(readlink -f "$script_merge_chrs" 2>/dev/null || realpath "$script_merge_chrs")"
  cd ${runroot}
  source config_shell.sh

  mkdir -p "\${LOG_PATH}/GWAS_Preprocessing" "${params.output_path}/GWAS_Preprocessing"

    bash "\$SCRIPT_PATH"

  # Stage the merged file back into runroot so NF can capture it
  ln -sf "${params.output_path}/GWAS_Preprocessing/${trait}_regenie_allchr.txt" "./${trait}_regenie_allchr.txt"
  """

  stub:
  """
  mkdir -p runroot
  : > runroot/${trait}_regenie_allchr.txt
  """
}
