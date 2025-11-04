/*
 * modules/druggability/non_disease_term_drug_profile.nf
 *
 * Purpose:
 *   Create non-disease-term drug profiles for gene-drug associations and
 *   publish results into the Druggability output folder.
 *
 * Inputs:
 *   - tuple(runroot_path, prioritised_genes_csv)
 *
 * Outputs:
 *   - nondiseaseterm_gene-drug_OT_FDA_adverse_effects.csv
 */
nextflow.enable.dsl=2

process NON_DISEASE_TERM_DRUG_PROFILE {
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
    if [ -f "${params.output_path}/Druggability/nondiseaseterm_gene-drug_OT_FDA_adverse_effects.csv" ]; then
        echo "✓ Non-disease term drug profiles already exist, skipping"
        exit 0
    fi

    # Change to runroot directory where config files are located
    cd "${runroot_path}"

    echo "Running non-disease term drug profile generation..."
    Rscript "${projectDir}/bin/druggability/non_disease_term_drug_profile.R"
    echo "✓ Non-disease term drug profiles complete"
    """

    stub:
    """
    echo "Stub: NON_DISEASE_TERM_DRUG_PROFILE"
    """
}
