/*
 * modules/prioritisation/gene_prioritisation.nf
 *
 * Purpose:
 *   Wrap the R-based gene prioritisation step. Stages optional inputs
 *   (PoPS, custom coloc) into the process workdir and invokes
 *   `bin/genes/gene_prioritisation.R` which writes Prioritised_genes.csv.
 *
 * Notes:
 *   - The custom-coloc CSV is accepted as a string (may point to a missing
 *     file); the process will create a harmless placeholder so the R script
 *     can run and set coloc-related flags to 'No' when appropriate.
 */
nextflow.enable.dsl = 2

process GENE_PRIORITISATION {
  label 'heavy'
  tag "all_traits"

  publishDir "${params.output_path}/Prioritisation",
             mode: 'copy',
             overwrite: true,
             pattern: "Prioritised_genes.csv"

  input:
  // Consumes: (runroot_path, PoPS top genes *path as STRING*, CUSTOM coloc csv as STRING)
  // We accept the coloc CSV as a value (string) rather than a 'path' so the file is optional
  // and may be missing; the process will create a harmless placeholder file when absent.
  tuple val(runroot_path), val(pops_top_genes_path), val(custom_coloc_results)

  output:
  // Emits: (runroot_path, Prioritised_genes.csv) from the task work dir
  tuple val(runroot_path), path("Prioritised_genes.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/genes/gene_prioritisation.R"
  // Use provided custom_trait_name when available; otherwise fall back to a safe
  // placeholder name so prioritisation can run when custom coloc is disabled.
  def custom_name = (params.containsKey('custom_trait_name') && params.custom_trait_name?.toString()?.trim()) ? params.custom_trait_name?.toString() : 'CUSTOM_DISABLED'
  def expected_coloc_basename = "${custom_name}_coloc_all_variants_pp4.csv"

  """
  set -euo pipefail

  WORK_DIR="\$PWD"
  SCRIPT_PATH="${script_path}"
  RUNROOT_PATH="${runroot_path}"
  REQUIRE_POPS="${params.containsKey('run_pops') ? (params.run_pops as boolean) : true}"
  POPS_TOP_PATH="${pops_top_genes_path}"

  # Custom coloc CSV is staged by Nextflow in the TASK WORK DIR
  COLOC_SRC_PATH="${custom_coloc_results}"

  # 1) Validate staged custom-coloc file presence (may be empty; R handles it).
  # If the coloc CSV is missing, create a harmless empty placeholder so the R script can run.
  if [ ! -e "\${COLOC_SRC_PATH}" ]; then
    ALT_PATH="${params.output_path}/Colocalisation/${expected_coloc_basename}"
    if [ -e "\${ALT_PATH}" ]; then
      # prefer the canonical results location if present
      COLOC_SRC_PATH="\${ALT_PATH}"
    else
      echo "WARN: CUSTOM coloc file missing at staged path: \${COLOC_SRC_PATH}" >&2
      echo "      and not found at: ${params.output_path}/Colocalisation/${expected_coloc_basename}" >&2
      echo "      Creating empty placeholder so prioritisation can continue." >&2
  # create a minimal placeholder file in the work dir with the expected basename
  # include a minimal header (Locus_name) so downstream fread() sees a well-formed CSV
  printf "Locus_name\n" > "\${WORK_DIR}/${expected_coloc_basename}"
      COLOC_SRC_PATH="${expected_coloc_basename}"
    fi
  fi

  # 2) Validate PoPS presence if required (must be non-empty when run_pops=true)
  if [ "\${REQUIRE_POPS}" = "true" ]; then
    # Accept the run-level placeholder `/dev/null` (used when PoPS output is
    # missing) or an empty file. In those cases warn and proceed without PoPS
    # rather than failing the whole pipeline. This makes the pipeline robust to
    # timing/publish issues where PoPS output may not be present yet.
    if [ "\${POPS_TOP_PATH}" = "/dev/null" ] || [ ! -s "\${POPS_TOP_PATH}" ]; then
      echo "WARN: PoPS required but missing/empty at: \${POPS_TOP_PATH}. Proceeding without PoPS." >&2
      # Flip the flag so downstream R logic knows PoPS is not available
      REQUIRE_POPS="false"
      POPS_TOP_PATH="/dev/null"
    fi
  else
    echo "INFO: run_pops=false — proceeding without PoPS file."
  fi

  # 3) Place/symlink the custom-coloc file INSIDE runroot BEFORE cd, using absolute paths
  DEST_COLOC_PATH="\${RUNROOT_PATH}/${expected_coloc_basename}"

  # Remove any existing file or symlink at the destination so we can create
  # a fresh link or copy. This avoids errors such as
  # "cp: not writing through dangling symlink ..." when a stale symlink exists.
  rm -f "\${DEST_COLOC_PATH}" || true

    # Prefer a lightweight symlink; fall back to copy if symlink fails (e.g., FS restrictions)
    if [ ! -e "\${DEST_COLOC_PATH}" ]; then
      if ln -sfn "\${WORK_DIR}/\${COLOC_SRC_PATH}" "\${DEST_COLOC_PATH}" 2>/dev/null; then
        :
      else
        # Safer copy strategy: copy to a temp file inside the runroot directory
        # then atomically move into place. This replaces any existing dangling
        # symlink at the destination rather than attempting to write through it.
  TMP_DEST="\$(dirname \"\${DEST_COLOC_PATH}\")/.tmp_coloc_\$\$"
  mkdir -p "\$(dirname \"\${DEST_COLOC_PATH}\")"
        if cp -f "\${WORK_DIR}/\${COLOC_SRC_PATH}" "\${TMP_DEST}" 2>/dev/null || \
           cp -f "\${COLOC_SRC_PATH}" "\${TMP_DEST}" 2>/dev/null; then
          mv -f "\${TMP_DEST}" "\${DEST_COLOC_PATH}"
        else
          # As a last resort write a minimal placeholder and move it into place
          printf "Locus_name\n" > "\${TMP_DEST}"
          mv -f "\${TMP_DEST}" "\${DEST_COLOC_PATH}"
        fi
      fi
    fi

  # Robustness: if the above link/copy didn't materialise the file in runroot,
  # attempt an explicit copy from the work dir or create a minimal file directly
  # inside runroot so downstream R sees a file for fread().
    if [ ! -e "\${DEST_COLOC_PATH}" ]; then
      if [ -e "\${WORK_DIR}/\${COLOC_SRC_PATH}" ]; then
        TMP_DEST="\$(dirname \"\${DEST_COLOC_PATH}\")/.tmp_coloc_\$\$"
        mkdir -p "\$(dirname \"\${DEST_COLOC_PATH}\")"
        if cp -f "\${WORK_DIR}/\${COLOC_SRC_PATH}" "\${TMP_DEST}" 2>/dev/null; then
          mv -f "\${TMP_DEST}" "\${DEST_COLOC_PATH}"
        fi
    elif [ -e "\${COLOC_SRC_PATH}" ]; then
        TMP_DEST="\$(dirname \"\${DEST_COLOC_PATH}\")/.tmp_coloc_\$\$"
        mkdir -p "\$(dirname \"\${DEST_COLOC_PATH}\")"
        if cp -f "\${COLOC_SRC_PATH}" "\${TMP_DEST}" 2>/dev/null; then
          mv -f "\${TMP_DEST}" "\${DEST_COLOC_PATH}"
        fi
    else
      # As a last resort create a minimal placeholder directly in runroot
      mkdir -p "\$(dirname \"\${DEST_COLOC_PATH}\")"
      printf "Locus_name\n" > "\${DEST_COLOC_PATH}"
    fi
  fi

  # Legacy CAD filename compatibility (only create if different and missing)
  if [ "${expected_coloc_basename}" != "CAD_coloc_all_variants_pp4.csv" ]; then
    LEGACY_CAD_PATH="\${RUNROOT_PATH}/CAD_coloc_all_variants_pp4.csv"
    if [ ! -e "\${LEGACY_CAD_PATH}" ]; then
      if ln -sfn "\${WORK_DIR}/\${COLOC_SRC_PATH}" "\${LEGACY_CAD_PATH}" 2>/dev/null; then
        :
      else
          TMP_LEG="\$(dirname \"\${LEGACY_CAD_PATH}\")/.tmp_coloc_\$\$"
          mkdir -p "\$(dirname \"\${LEGACY_CAD_PATH}\")"
          if cp -f "\${COLOC_SRC_PATH}" "\${TMP_LEG}" 2>/dev/null; then
            mv -f "\${TMP_LEG}" "\${LEGACY_CAD_PATH}"
          fi
      fi
    fi
  fi

  # 4) Now enter runroot
  cd "\${RUNROOT_PATH}"

  # Sanity check that the expected file is visible from here (existence only; may be empty)
  if [ ! -e "${expected_coloc_basename}" ]; then
    echo "ERROR: Expected coloc file not found inside runroot: \${PWD}/${expected_coloc_basename}" >&2
    echo "       It should have been linked/copied from: \${WORK_DIR}/\${COLOC_SRC_PATH}" >&2
    exit 1
  fi

  # 5) Ensure target output folder exists
  mkdir -p "${params.output_path}/Prioritisation"

  # 6) Run the prioritisation R script (handles empty coloc files gracefully)
  Rscript "\${SCRIPT_PATH}"

  OUT_SRC="${params.output_path}/Prioritisation/Prioritised_genes.csv"

  # 7) Emit output
  if [ -s "\${OUT_SRC}" ]; then
    cp -f "\${OUT_SRC}" "\${WORK_DIR}/Prioritised_genes.csv" || ln -sf "\${OUT_SRC}" "\${WORK_DIR}/Prioritised_genes.csv"
  else
    echo "ERROR: Expected output not produced: \${OUT_SRC}" >&2
    exit 2
  fi

  ls -l "\${WORK_DIR}/Prioritised_genes.csv"
  """
}
