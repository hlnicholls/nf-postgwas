/*
 * modules/druggability/gene_drug_target_profiles.nf
 *
 * Purpose:
 *   Produce gene-drug target profiles for prioritised genes and publish the
 *   aggregated CSV into the Druggability output folder.
 *
 * Inputs:
 *   - tuple(runroot_path, prioritised_genes_csv)
 *
 * Outputs:
 *   - Druggable_gene_targets.csv
 */
nextflow.enable.dsl=2

process GENE_DRUG_TARGET_PROFILES {
    label 'light'
    tag "all_traits"
    publishDir "${params.output_path}/Druggability", mode: 'copy', overwrite: true

    input:
    tuple val(runroot_path), path(prioritised_genes_csv)

    output:
    tuple val(runroot_path), path(prioritised_genes_csv)

    script:
    """
    set -euo pipefail

    # Skip if output already exists
    if [ -f "${params.output_path}/Druggability/Druggable_gene_targets.csv" ]; then
        echo "✓ Gene drug target profiles already exist, skipping"
        exit 0
    fi

    # Change to runroot directory where config files are located
    cd "${runroot_path}"

    echo "Running gene-drug target profile generation..."
    Rscript "${projectDir}/bin/druggability/gene_drug_target_profiles.R"
    echo "✓ Gene-drug target profiles complete"
    """

    stub:
    """
    echo "Stub: GENE_DRUG_TARGET_PROFILES"
    """
}
