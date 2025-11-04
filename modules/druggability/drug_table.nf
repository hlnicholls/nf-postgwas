/*
 * modules/druggability/drug_table.nf
 *
 * Purpose:
 *   Generate a druggability table from the prioritised gene list and publish
 *   the resulting CSV into the Druggability output folder.
 *
 * Inputs:
 *   - tuple(runroot_path, prioritised_genes_csv)
 *
 * Outputs:
 *   - Druggability_table.csv (published to ${params.output_path}/Druggability)
 */
nextflow.enable.dsl=2

process DRUG_TABLE {
    label 'medium'
    tag "all_traits"
    publishDir "${params.output_path}/Druggability", mode: 'copy', overwrite: true

    input:
    // We cd to this path (runroot)
    tuple val(runroot_path), path(prioritised_genes_csv)

    output:
    // Pass the same tuple forward for downstream modules
    tuple val(runroot_path), path(prioritised_genes_csv)

    script:
    """
    set -euo pipefail

    # Skip if output already exists
    if [ -f "${params.output_path}/Druggability/Druggability_table.csv" ]; then
        echo "✓ Drug table already exists, skipping"
        exit 0
    fi

    # Ensure output directory exists (publishDir will copy CSVs from work dir)
    mkdir -p "${params.output_path}/Druggability"

    # Change to runroot directory where config files are located
    cd "${runroot_path}"

    echo "Running drug table generation..."
    Rscript "${projectDir}/bin/druggability/drug_table.R"
    echo "✓ Drug table generation complete"
    """

    stub:
    """
    echo "Stub: DRUG_TABLE"
    """
}
