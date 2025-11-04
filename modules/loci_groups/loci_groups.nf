/*
 * modules/loci_groups/loci_groups.nf
 *
 * Purpose:
 *   Group nearby loci into blocks using distance and LD thresholds and
 *   produce an 'All_loci_blocks.csv' file for downstream PoPS and plotting.
 *
 * Inputs:
 *   - tuple(runroot_path, loci_with_ld_genes)
 *
 * Outputs:
 *   - All_loci_blocks.csv (published to Loci_Preprocessing)
 */
nextflow.enable.dsl=2

process LOCI_GROUPS {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Loci_Preprocessing",
             mode: 'copy',
             overwrite: true,
             pattern: "All_loci_blocks.csv"

  input:
  tuple val(runroot_path), path(loci_with_ld_genes)

  output:
  tuple val(runroot_path), path("All_loci_blocks.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/loci/loci_groups_1mb_4mbLD.R"
  """
  set -euo pipefail

  # Stage inputs so the R script sees the expected filename in runroot
  test -s "all_traits_loci_38_with_ld_genes.txt" >/dev/null
  cp all_traits_loci_38_with_ld_genes.txt "${runroot_path}/all_traits_loci_38_with_ld_genes.txt"

  WORK_DIR=\$PWD
  cd "${runroot_path}"

  # Ensure config exists
  [ -f config_R.R ] || { echo "Missing config_R.R in ${runroot_path}" >&2; exit 1; }

  # --- Sanitize config_R.R -------------------------------------------------
  # Remove any absolute (unquoted) reassignment lines injected by other steps
  # e.g. `var_file <- /abs/path/...`
  sed -i '/^[[:space:]]*var_file[[:space:]]*<-[[:space:]]*\\//d' config_R.R

  # Re-assert portable var_file definition (based on loci_output_path)
  # This stays portable across machines/paths.
  awk '
    BEGIN{printed=0}
    {print}
    END{
      if (!printed) {
        print "var_file <- file.path(loci_output_path, \\"all_traits_loci_38_with_ld_genes.txt\\")"
      }
    }
  ' config_R.R > config_R.R.tmp && mv config_R.R.tmp config_R.R
  # -------------------------------------------------------------------------

  # Run the loci grouping script (it sources config_R.R)
  Rscript "${script_path}"

  # Copy output back to work directory for Nextflow to publish
  cp "${params.output_path}/Loci_Preprocessing/All_loci_blocks.csv" "\$WORK_DIR/"
  """
}
