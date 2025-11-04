/*
 * modules/coloc/custom_coloc.nf
 *
 * Purpose:
 *   Run a user-provided custom colocalisation script and process its outputs
 *   into a single CSV file consumable by downstream prioritisation steps.
 *
 * Inputs:
 *   - tuple runroot_path, all_loci, script_custom_coloc, script_process_custom
 *
 * Outputs:
 *   - ${params.custom_trait_name}_coloc_all_variants_pp4.csv (placeholder created if missing)
 */
nextflow.enable.dsl=2

process CUSTOM_COLOC {
  label 'heavy'
  tag "all_traits"

  publishDir "${params.output_path}/Colocalisation",
             mode: 'copy',
             overwrite: true,
             pattern: "${params.custom_trait_name}_coloc_*.csv"

  input:
  tuple val(runroot_path), path(all_loci), path(script_custom_coloc), path(script_process_custom)

  output:
  // keep this non-optional so downstream gets a token emission
  path("${params.custom_trait_name}_coloc_all_variants_pp4.csv")

  script:
  """
  set -euo pipefail
  set -x

  WORK_DIR="\$PWD"
  echo "Staged files in \$WORK_DIR:"
  ls -l

  SCRIPT_CUSTOM_PATH="${script_custom_coloc}"
  SCRIPT_PROCESS_PATH="${script_process_custom}"
  ALL_LOCI_PATH="${all_loci}"

  [ -s "\$SCRIPT_CUSTOM_PATH" ]  || { echo "ERROR: custom script missing: \$SCRIPT_CUSTOM_PATH"; exit 2; }
  [ -s "\$SCRIPT_PROCESS_PATH" ] || { echo "ERROR: process script missing: \$SCRIPT_PROCESS_PATH"; exit 2; }
  [ -s "\$ALL_LOCI_PATH" ]       || { echo "ERROR: loci CSV missing: \$ALL_LOCI_PATH"; exit 2; }

  COLOC_DIR="${params.output_path}/Colocalisation"
  mkdir -p "\$COLOC_DIR"

  # Count phenotypes present in the loci CSV (if any)
  EXPECTED_PHENOS=\$(awk -F',' 'NR>1 {print \$1}' "\$ALL_LOCI_PATH" | sort -u | wc -l | tr -d ' ')
  echo "EXPECTED_PHENOS=\$EXPECTED_PHENOS"

  # Pre-count existing RDS results (excluding *noPP4*)
  RDS_FILES=\$(find "\$COLOC_DIR" -maxdepth 1 -type f -name "coloc_${params.custom_trait_name}_res_all_*.rds" ! -name "*noPP4*" | wc -l | tr -d ' ')
  echo "Found \$RDS_FILES existing RDS files in \$COLOC_DIR"

  # Work from runroot so config_R.R & relative paths resolve
  cd "${runroot_path}"
  ln -sf "\$WORK_DIR/\$SCRIPT_CUSTOM_PATH"   ./custom_coloc.R
  ln -sf "\$WORK_DIR/\$SCRIPT_CUSTOM_PATH"   ./Custom_coloc.R
  ln -sf "\$WORK_DIR/\$SCRIPT_PROCESS_PATH"  ./process_custom_coloc_output.R
  ln -sf "\$WORK_DIR/\$SCRIPT_PROCESS_PATH"  ./process_Custom_coloc_output.R
  ln -sf "\$WORK_DIR/\$ALL_LOCI_PATH"        ./All_loci_ungrouped.csv

  # Run only if we don't already have >= expected RDS (and at least one expected)
  if [ "\$EXPECTED_PHENOS" -gt 0 ] && [ "\$RDS_FILES" -lt "\$EXPECTED_PHENOS" ]; then
    echo "Running custom coloc..."
    if [ -s ./custom_coloc.R ]; then
      Rscript --vanilla ./custom_coloc.R 2>&1 | tee "\$WORK_DIR/custom_coloc.stdout.log"
    else
      Rscript --vanilla ./Custom_coloc.R 2>&1 | tee "\$WORK_DIR/custom_coloc.stdout.log"
    fi
  else
    echo "Skipping custom coloc run (nothing expected or already satisfied)"
  fi

  echo "Processing custom coloc outputs..."
  if [ -s ./process_custom_coloc_output.R ]; then
    Rscript --vanilla ./process_custom_coloc_output.R 2>&1 | tee "\$WORK_DIR/custom_coloc.process.stdout.log"
  else
    Rscript --vanilla ./process_Custom_coloc_output.R 2>&1 | tee "\$WORK_DIR/custom_coloc.process.stdout.log"
  fi

  # Ensure the declared Nextflow output file is PRESENT in the work dir
  OUT_BASENAME="${params.custom_trait_name}_coloc_all_variants_pp4.csv"
  SRC_PATH="\${COLOC_DIR}/\${OUT_BASENAME}"

  # If the processor wrote it to COLOC_DIR, copy it back regardless of size
  if [ -f "\$SRC_PATH" ]; then
    cp -f "\$SRC_PATH" "\$WORK_DIR/\$OUT_BASENAME" || true
  fi

  # If it still doesn't exist in the work dir (e.g. processor didn’t create it),
  # create an empty placeholder so Nextflow can collect & publish it.
  if [ ! -f "\$WORK_DIR/\$OUT_BASENAME" ]; then
    : > "\$WORK_DIR/\$OUT_BASENAME"
  fi
  """
}
