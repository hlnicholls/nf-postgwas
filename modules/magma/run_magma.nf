/*
 * modules/magma/run_magma.nf
 *
 * Purpose:
 *   Run MAGMA gene-based analysis on prepped per-trait association files.
 *
 * Inputs:
 *   - tuple runroot_path, trait, _is_mtag, prepped_file
 *
 * Outputs:
 *   - magma_${trait}.genes.out (published under MAGMA/Results)
 */
nextflow.enable.dsl=2

process RUN_MAGMA {
  label 'magma'
  tag "${trait}"

  publishDir "${params.output_path}/MAGMA/Results",
             mode: 'copy',
             overwrite: true,
             pattern: "magma_*.genes.out"

  input:
  tuple val(runroot_path), val(trait), val(_is_mtag), path(prepped_file)

  output:
  tuple val(trait), path("magma_*.genes.out")

  script:
  // Config values with defaults
  // Keep data assets in the external databases path
  def magma_bfile = "${params.databases}/MAGMA/1000G.EUR"
  def magma_annot = "${params.databases}/MAGMA/magma_0kb.genes.annot"

  // Use the in-image MAGMA binary on PATH; allow override via --magma_bin
  def magma_bin   = params.magma_bin ?: "magma"

  def output_prefix = "magma_${trait}"
  // Use the staged file from work directory
  def input_file = "\$WORK_DIR/${prepped_file}"

  """
  set -euo pipefail

  WORK_DIR=\$PWD
  cd ${runroot_path}

  # Ensure results dir exists
  mkdir -p ${params.output_path}/MAGMA/Results

  # Verify MAGMA binary is available
  if ! command -v ${magma_bin} >/dev/null 2>&1; then
    echo "ERROR: MAGMA binary '${magma_bin}' not found on PATH. Ensure it is installed in the image." >&2
    exit 127
  fi

  # --- Run MAGMA gene-based analysis ---
  ${magma_bin} \\
    --bfile ${magma_bfile} \\
    --gene-annot ${magma_annot} \\
    --pval ${input_file} ncol=N \\
    --gene-model snp-wise=mean \\
    --out ${params.output_path}/MAGMA/Results/${output_prefix}

  # Stage outputs back to work dir
  cp ${params.output_path}/MAGMA/Results/${output_prefix}.genes.out "\$WORK_DIR/" || true
  """
}
