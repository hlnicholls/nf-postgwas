/*
 * modules/enrichment/enrichr_annotation.nf
 *
 * Purpose:
 *   Run the EnrichR gene annotation step on genes output from loci processing
 *   and publish annotated gene tables to the Enrichment output folder.
 *
 * Inputs:
 *   - tuple runroot_path, loci_with_ld_genes, script_enrichr
 *
 * Outputs:
 *   - All_genes_annotated_with_EnrichR_*.csv
 */
nextflow.enable.dsl=2

process ENRICHR_ANNOTATION {
  label 'heavy'
  tag "all_traits"

  publishDir "${params.output_path}/Enrichment",
             mode: 'copy',
             overwrite: true,
             pattern: "All_genes_annotated_with_EnrichR_*.csv"

  input:
  tuple val(runroot_path), path(loci_with_ld_genes), path(script_enrichr)

  output:
  path "All_genes_annotated_with_EnrichR_*.csv"

  script:
  """
  set -euo pipefail
  
  # Store the work directory
  WORK_DIR=\$PWD
  
  # Resolve script path before changing directories
  SCRIPT_PATH="\$(readlink -f "$script_enrichr" 2>/dev/null || realpath "$script_enrichr")"
  
  # Run the script from within the runroot
  cd ${runroot_path}
  
  # Ensure enrichment directory exists
  mkdir -p "${params.output_path}/Enrichment"
  
  # Create the expected directory structure for here::i_am()
  mkdir -p Enrichment

  # Create a symlink to the script in the expected location
  ln -sf "\$SCRIPT_PATH" Enrichment/Step2_Enrichr_Gene_annotation.R
  
  # Run the enrichr annotation script — let non-zero exit propagate and fail the process
  Rscript "\$SCRIPT_PATH"

  # Copy outputs back to work directory (will error if outputs are missing)
  cp "${params.output_path}/Enrichment"/All_genes_annotated_with_EnrichR_*.csv "\$WORK_DIR/"
  """
}
