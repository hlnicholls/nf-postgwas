/*
 * modules/prep/link_premade_gwas.nf
 *
 * Purpose:
 *   If the user supplied a pre-made merged GWAS file, copy/atomically-move
 *   it into the canonical pipeline results directory and mirror into runroot
 *   so subsequent PREP steps can proceed unchanged.
 *
 * Inputs:
 *   - tuple(trait, runroot, merged_file)
 *
 * Outputs:
 *   - <trait>_regenie_allchr.txt (mirrored into runroot)
 */
nextflow.enable.dsl=2

process LINK_PREMADE_GWAS {
  tag { trait }
  label 'light'

  input:
  tuple val(trait), path(runroot), path(merged_file)

  output:
  tuple val(trait), path("runroot/${trait}_regenie_allchr.txt")

  /*
   * NOTE:
   * - Use Groovy interpolation for Nextflow params: ${params.output_path}
   * - Escape only real shell vars ($TMP, $DEST_RESULTS, $$, etc.)
   */
  script:
  """
  set -euo pipefail

  # Canonical results location (same as MERGE_CHRS writes to)
  RES_DIR="${params.output_path}/GWAS_Preprocessing"
  DEST_RESULTS="\$RES_DIR/${trait}_regenie_allchr.txt"
  DEST_RUNROOT="runroot/${trait}_regenie_allchr.txt"

  # Ensure directories exist (don't clobber staged runroot)
  mkdir -p "\$RES_DIR"
  [ -d "runroot" ] || mkdir -p "runroot"

  # Resolve absolute source path robustly (Linux/Mac)
  SRC="\$(readlink -f "${merged_file}" 2>/dev/null || realpath "${merged_file}" 2>/dev/null || echo "${merged_file}")"

  # Copy into canonical results folder (atomic-ish to avoid partials)
  TMP="\${DEST_RESULTS}.tmp.\$\$"
  cp -f "\$SRC" "\$TMP"
  mv -f "\$TMP" "\$DEST_RESULTS"

  # Mirror into runroot for NF capture; prefer symlink, else copy
  if ln -sfn "\$DEST_RESULTS" "\$DEST_RUNROOT" 2>/dev/null; then
    :
  else
    cp -f "\$DEST_RESULTS" "\$DEST_RUNROOT"
  fi
  """

  stub:
  """
  mkdir -p "runroot" "${params.output_path}/GWAS_Preprocessing"
  : > "${params.output_path}/GWAS_Preprocessing/${trait}_regenie_allchr.txt"
  ln -sfn "${params.output_path}/GWAS_Preprocessing/${trait}_regenie_allchr.txt" "runroot/${trait}_regenie_allchr.txt" || \
    : > "runroot/${trait}_regenie_allchr.txt"
  """
}
