/*
 * modules/plots/regional_plots.nf
 *
 * Purpose:
 *   Produce regional locus plots (locuszoom-style) for loci using LD
 *   annotations. Writes PNGs to the GWAS_Plots/Regional_plots directory
 *   and emits a 'regional_plots_done.txt' sentinel on completion.
 *
 * Inputs:
 *   - tuple (runroot_path, ld_annotated, script_locuszoom, script_makeRegionalPlot)
 *
 * Outputs:
 *   - regional_plots_done.txt
 */
nextflow.enable.dsl=2

process REGIONAL_PLOTS {
  label 'medium'
  tag 'all_traits'
  errorStrategy 'ignore'

  publishDir "${params.output_path}/GWAS_Plots/Regional_plots",
             mode: 'copy',
             overwrite: true,
             pattern: "*.png"

  input:
  tuple val(runroot_path), path(ld_annotated), path(script_locuszoom), path(script_makeRegionalPlot)

  output:
  path("regional_plots_done.txt")

  script:
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  # Resolve script paths
  SCRIPT_LZ=\$(readlink -f "${script_locuszoom}" 2>/dev/null || realpath "${script_locuszoom}")
  SCRIPT_MRP=\$(readlink -f "${script_makeRegionalPlot}" 2>/dev/null || realpath "${script_makeRegionalPlot}")

  RUNROOT_PATH="${runroot_path}"
  OUTDIR="${params.output_path}/GWAS_Plots/Regional_plots"

  [ -d "\$RUNROOT_PATH" ] || { echo "Missing runroot: \$RUNROOT_PATH" >&2; exit 1; }
  [ -s "\$SCRIPT_LZ" ]    || { echo "Missing plot script: \$SCRIPT_LZ" >&2; exit 1; }
  [ -s "\$SCRIPT_MRP" ]   || { echo "Missing makeRegionalPlot.R: \$SCRIPT_MRP" >&2; exit 1; }

  # Stage scripts; do not touch config_R.R content
  cp "\$SCRIPT_LZ"  "\$RUNROOT_PATH/locuszoom_plots.R"
  cp "\$SCRIPT_MRP" "\$RUNROOT_PATH/makeRegionalPlot.R"

  cd "\$RUNROOT_PATH"
  [ -f config_R.R ] || { echo "Missing config_R.R in \$RUNROOT_PATH" >&2; exit 1; }

  # --- Sanitize config_R.R -------------------------------------------------
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

  # Ensure output directory exists
  mkdir -p "\$OUTDIR"

  # Run the regional plotter (sources config_R.R)
  Rscript "locuszoom_plots.R"

  echo "Regional plots completed" > "\$WORK_DIR/regional_plots_done.txt"
  """
}
