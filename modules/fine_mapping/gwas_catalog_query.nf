/*
 * modules/fine_mapping/gwas_catalog_query.nf
 *
 * Purpose:
 *   Query the GWAS Catalog for credible-set variant overlaps and publish
 *   pleiotropy summary CSVs into the Pleiotropy output folder.
 *
 * Inputs:
 *   - tuple(runroot_path, annotated_credsets)
 *
 * Outputs:
 *   - GWAS catalog query CSV(s) under ${params.output_path}/Pleiotropy
 */
nextflow.enable.dsl=2

process GWAS_CATALOG_QUERY {
  label 'heavy'
  tag "all_traits"

  publishDir "${params.output_path}/Pleiotropy",
             mode: 'copy',
             overwrite: true,
             pattern: "*.csv"

  input:
  tuple val(runroot_path), path(annotated_credsets)

  output:
  tuple val(runroot_path), path("*.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/credsets/GWAS_Catalog_Query.R"
  def outputDir = file("${params.output_path}/Pleiotropy")
  
  // Check if GWAS Catalog query results already exist
  def resultsExist = outputDir.exists() && 
                     outputDir.list().any { it.contains('GWAS_catalog_query') }
  
  if (resultsExist) {
    """
    # GWAS Catalog query already exists - copy to work directory
  cp "${params.output_path}/Pleiotropy"/GWAS_catalog_query*.csv . 2>/dev/null || true
    """
  } else {
    """
    set -euo pipefail

    WORK_DIR=\$PWD

    # Touch input to ensure it exists
    test -s "${annotated_credsets}" > /dev/null

    cd ${runroot_path}

    # Ensure output directory exists
  mkdir -p "${params.output_path}/Pleiotropy"

    # Use shared config (created by CREATE_CONFIG_SHIMS)
    Rscript ${script_path}

    # Stage outputs back to work dir
  cp "${params.output_path}/Pleiotropy"/*.csv "\$WORK_DIR/" 2>/dev/null || true
    """
  }
}
