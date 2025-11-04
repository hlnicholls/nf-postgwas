/*
 * modules/ld/calc_ld.nf
 *
 * Purpose:
 *   Calculate per-trait LD matrices using the provided LD helper script.
 *   This process runs inside each trait's runroot and emits a small
 *   sentinel file when the per-trait LD computation completes.
 *
 * Inputs:
 *   - tuple(trait, runroot_dir, rsids_txt, loci_done_flag, ld_script)
 *
 * Outputs:
 *   - <runroot>/<trait>_ld_raw.txt (sentinel)
 */
nextflow.enable.dsl=2

process CALC_LD {
  label 'heavy'
  tag "${trait}"

  input:
  tuple val(trait),
        path(runroot_dir),
        path(rsids_txt),
        path(loci_done_flag),
        path(ld_script)

  output:
  tuple val(trait),
        path("${runroot_dir}/${trait}_ld_raw.txt")

  // nothing to publish here; *.ld are written to global OUTPUT_PATH by your bash

  script:
  // Use WORK_DIR to reference staged inputs from *inside* runroot
  def ld_script_path = "\$WORK_DIR/${ld_script}"

  """
  set -euo pipefail

  WORK_DIR=\$PWD
  RUNROOT="\$PWD/${runroot_dir}"

  # Optional: let users override PLINK memory via Nextflow param (falls back inside calculate_LD.sh)
  ${ params.plink_mem_mb ? "export PLINK_MEM_MB='${params.plink_mem_mb}'" : ":" }

  # sanity: config from shims must be present
  if [ ! -f "\$RUNROOT/config_shell.sh" ]; then
    echo "config_shell.sh missing in \$RUNROOT" >&2
    exit 1
  fi

  # tiny debug (kept quiet unless there’s a failure)
  echo "CWD: \$PWD" >&2
  echo "RUNROOT: \$RUNROOT" >&2
  echo "LD script (staged): ${ld_script_path}" >&2
  ${ params.plink_mem_mb ? "echo \"PLINK_MEM_MB override: \$PLINK_MEM_MB\" >&2" : ":" }
  ls -l "\$RUNROOT" >/dev/null 2>&1 || true
  ls -l "${ld_script_path}" >/dev/null 2>&1 || true

  # Run from inside runroot so ./config_shell.sh is in scope for your bash
  (
    cd "\$RUNROOT"
    bash "${ld_script_path}"
  )

  # Emit a sentinel so Nextflow can track completion
  : > "\$RUNROOT/${trait}_ld_raw.txt"
  """
}
