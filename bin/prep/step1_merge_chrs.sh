# 1_Loci_Preprocessing/1_QC_and_IDs/Step0_merge_chrs.sh (make idempotent + safer)
#!/usr/bin/env bash
set -euo pipefail

# When called from Nextflow, config_shell.sh should be in the current directory
# Use ./ prefix to ensure we source from current working directory, not from script's directory
if [[ -f "./config_shell.sh" ]]; then
  source "./config_shell.sh"
else
  SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  BASE_DIR="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
  CONFIG_FILE="${BASE_DIR}/config_shell.sh"
  source "${CONFIG_FILE}"
fi

mkdir -p "${LOG_PATH}/GWAS_Preprocessing" "${COMBINED_REGENIE_PATH}"

log_file="${LOG_PATH}/GWAS_Preprocessing/1_merge_regenie_file_chrs.log"

{
  cd "${REGENIE_FOLDER}"

  for trait in "${TRAITS[@]}"; do
    out_file="${COMBINED_REGENIE_PATH}/${trait}_regenie_allchr.txt"

    # Skip if already merged and non-empty (idempotent)
    if [[ -s "${out_file}" ]]; then
      echo "Skipping ${trait} (exists: ${out_file})"
      continue
    fi

    echo "Processing phenotype: ${trait}"

    # Find a sample file for header
    sample_file="$(ls | grep -E "(^|_)${trait}(_|\.regenie$)" | grep -v "${trait}_" | head -n 1 || true)"
    if [[ -z "${sample_file}" ]]; then
      echo "No files found for phenotype: ${trait}"
      continue
    fi

    # Write header (normalize spaces->tabs)
    head -n 1 "${sample_file}" | tr -s ' ' '\t' > "${out_file}"

    echo "Combining chromosome files for ${trait}:"
    # Pick all chr files for this trait (excluding already-merged names)
    while IFS= read -r f; do
      echo "  Adding file: ${f}"
      tail -n +2 "${f}" | tr -s ' ' '\t' >> "${out_file}"
    done < <(ls | grep -E "(^|_)${trait}(_|\.regenie$)" | grep -v "${trait}_assoc_regenie_allchr\.txt")

    echo "Completed combining files for ${trait}. Output: ${out_file}"
  done
} | tee -a "${log_file}"
