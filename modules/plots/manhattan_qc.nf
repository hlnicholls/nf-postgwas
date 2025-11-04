/*
 * modules/plots/manhattan_qc.nf
 *
 * Purpose:
 *   Generate per-trait Manhattan and QC plots using prepped GWAS inputs.
 *   The underlying R scripts write PNGs directly into the pipeline
 *   output directory; this process emits a small sentinel file per-trait
 *   so downstream steps can depend on completion without copying images.
 *
 * Inputs:
 *   - tuple (trait, runroot, rsids_file, script_qc_plots, script_manhattan)
 *
 * Outputs:
 *   - manhattan_qc_done_<trait>.txt
 */
nextflow.enable.dsl=2

process MANHATTAN_QC {
  label 'medium'
  tag { trait }

  /*
   * We intentionally do NOT publish *.png from the work dir.
   * The R scripts write directly into:
  *   ${params.output_path}/GWAS_Plots/QQ_plots
  *   ${params.output_path}/GWAS_Plots/QC_plots
   * We leave the files there and just emit a sentinel per trait.
   */

  input:
  tuple val(trait), path(runroot), path(rsids_file), path(script_qc_plots), path(script_manhattan)

  output:
  // Emit a small sentinel so downstream can depend on completion if needed
  tuple val(trait), path("manhattan_qc_done_${trait}.txt")

  script:
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  # Resolve script paths before changing directories
  SCRIPT_QC=\$(readlink -f "${script_qc_plots}" 2>/dev/null || realpath "${script_qc_plots}")
  SCRIPT_MANHATTAN=\$(readlink -f "${script_manhattan}" 2>/dev/null || realpath "${script_manhattan}")

  # Target output directories used by the R scripts
  QQ_DIR="${params.output_path}/GWAS_Plots/QQ_plots"
  QC_DIR="${params.output_path}/GWAS_Plots/QC_plots"

  # Ensure base plot subdirectories exist for scripts that don't create parents
  mkdir -p "\$QQ_DIR"
  mkdir -p "\$QC_DIR"

  # --- Skip if per-trait images already exist in BOTH locations ----------------
  # --- Skip if per-trait images already exist in QQ, QC AND Manhattan ---
  MAN_DIR="${params.output_path}/GWAS_Plots"
  if compgen -G "\$QQ_DIR/*${trait}*.png" > /dev/null && \
     compgen -G "\$QC_DIR/*${trait}*.png" > /dev/null && \
     compgen -G "\$MAN_DIR/*Manhattan*${trait}*.png" > /dev/null ; then
      echo "✓ Manhattan/QC outputs for ${trait} already exist — skipping"
      echo "done" > "\$WORK_DIR/manhattan_qc_done_${trait}.txt"
      exit 0
  fi
  # ---------------------------------------------------------------------------

  cd "${runroot}"

  # Require config
  [ -f config_R.R ] || { echo "Missing config_R.R in ${runroot}" >&2; exit 1; }

  # --- Sanitize config_R.R -------------------------------------------------
  # Remove any absolute var_file definition and ensure a default relative path is present
  sed -i '/^[[:space:]]*var_file[[:space:]]*<-[[:space:]]*\\//d' config_R.R
  awk '
    BEGIN{have=0}
    /^[[:space:]]*var_file[[:space:]]*<-/ { have=1 }
    {print}
    END{
      if (!have) {
        print "var_file <- file.path(loci_output_path, \\"all_traits_loci_38_with_ld_genes.txt\\")"
      }
    }
  ' config_R.R > config_R.R.tmp && mv config_R.R.tmp config_R.R
  # -------------------------------------------------------------------------

  # Run both scripts (they source config_R.R and write into subfolders)
  "\$SCRIPT_QC"
  "\$SCRIPT_MANHATTAN"

  # Create a sentinel so Nextflow can track completion without copying images
  echo "done" > "\$WORK_DIR/manhattan_qc_done_${trait}.txt"
  """
}
