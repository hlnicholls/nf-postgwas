/*
 * modules/variant_annot/to_vcf.nf
 *
 * Purpose:
 *   Convert summary statistics in LD to a VCF suitable for VEP annotation
 *   and publish the combined VEP input VCF file under Variant_annotation.
 *
 * Inputs:
 *   - tuple(runroot_path, ld_genes_file)
 *
 * Outputs:
 *   - all_traits_in_ld_VEP_input.vcf
 */
nextflow.enable.dsl=2

process VARIANT_ANNOT_TO_VCF {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Variant_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "all_traits_in_ld_VEP_input.vcf"

  input:
  tuple val(runroot_path), path(ld_genes_file)

  output:
  tuple val(runroot_path), path("all_traits_in_ld_VEP_input.vcf")

  script:
  def script_path = "${workflow.projectDir}/bin/variant_annot/convert_sumstats_to_VCF.R"
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  # Touch input to ensure it exists
  test -s "${ld_genes_file}" > /dev/null

  cd ${runroot_path}

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Variant_annotation"

  # Use shared config (created by CREATE_CONFIG_SHIMS)
  Rscript ${script_path}

  # Stage output back to work dir
  cp "${params.output_path}/Variant_annotation/all_traits_in_ld_VEP_input.vcf" "\$WORK_DIR/"
  """
}
