/*
 * modules/ldsc/clean_ldsc_logs.nf
 *
 * Purpose:
 *   Clean and reformat LDSC log outputs into tabular correlation tables and
 *   publish cleaned tables for downstream reporting.
 *
 * Inputs:
 *   - tuple(trait, runroot, gencorr_log, script_clean_ldsc_logs)
 *
 * Outputs:
 *   - per-trait *_corr_table.txt and all_pheno_corr_table.txt
 */
nextflow.enable.dsl=2

process CLEAN_LDSC_LOGS {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Loci_Preprocessing/LDSC_results",
             mode: 'copy',
             overwrite: true,
             pattern: "*.txt"

  input:
  tuple val(trait), path(runroot), path(gencorr_log), path(script_clean_ldsc_logs)

  output:
  tuple path("*_corr_table.txt"), path("all_pheno_corr_table.txt")

  script:
  """
set -euo pipefail

# Keep a reference to the current (work) directory
WORK_DIR="\$PWD"

# Resolve absolute path to the staged R script *before* changing directories
# Use readlink -f if available; otherwise fallback to realpath
if command -v readlink >/dev/null 2>&1; then
  SCRIPT_PATH="\$(readlink -f "$script_clean_ldsc_logs")"
else
  SCRIPT_PATH="\$(realpath "$script_clean_ldsc_logs")"
fi

# Sanity check to help debugging if staging ever breaks
[ -f "\$SCRIPT_PATH" ] || { echo "ERROR: R script not found at \$SCRIPT_PATH"; ls -la; exit 2; }

# Move into runroot so any relative paths inside the R script (e.g. config_R.R) work
cd "${runroot}"

# Append LDSC-specific paths to existing config_R.R instead of overwriting it
printf '%s\n' \
  'ldsc_results_path <- file.path(output_path, "Loci_Preprocessing/LDSC_results")' \
  'ldsc_log_path_for_table <- file.path(output_path, "Loci_Preprocessing")' \
  >> config_R.R

# Run the LDSC log cleaner via Rscript using the absolute path
Rscript "\$SCRIPT_PATH"

# Copy outputs back to the process work directory so Nextflow can capture declared outputs
cp ${params.output_path}/Loci_Preprocessing/LDSC_results/*_corr_table.txt "\$WORK_DIR/" || true
cp ${params.output_path}/Loci_Preprocessing/LDSC_results/all_pheno_corr_table.txt "\$WORK_DIR/" || true
  """
}
