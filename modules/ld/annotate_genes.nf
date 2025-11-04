/*
 * modules/ld/annotate_genes.nf
 *
 * Purpose:
 *   Annotate compiled loci with gene annotations and LD proxies using the
 *   provided Python/R helper. Publishes a combined loci-with-genes file.
 *
 * Inputs:
 *   - tuple(runroot_path, ld_collated, script_annotate_genes)
 *
 * Outputs:
 *   - all_traits_loci_38_with_ld_genes.txt (published to Loci_Preprocessing)
 */
nextflow.enable.dsl=2

process ANNOTATE_LD_GENES {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Loci_Preprocessing",
             mode: 'copy',
             overwrite: true,
             pattern: "all_traits_loci_38_with_ld_genes.txt"

  input:
  tuple val(runroot_path), path(ld_collated), path(script_annotate_genes)

  output:
  tuple val(runroot_path), path("all_traits_loci_38_with_ld_genes.txt")

  script:
  """
  set -euo pipefail
  
  # Store the work directory before changing
  WORK_DIR=\$PWD
  
  # Resolve script path before changing directories
  SCRIPT_PATH="\$(readlink -f "$script_annotate_genes" 2>/dev/null || realpath "$script_annotate_genes")"
  
  cd ${runroot_path}

    "\$SCRIPT_PATH"

  # Copy the output file to the work directory
  cp "${params.output_path}/Loci_Preprocessing/all_traits_loci_38_with_ld_genes.txt" "\$WORK_DIR/all_traits_loci_38_with_ld_genes.txt"
  """
}
