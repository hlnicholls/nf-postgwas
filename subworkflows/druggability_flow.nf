/*
 * subworkflows/druggability_flow.nf
 *
 * Purpose:
 *   Run a small chain of druggability analyses (drug table, gene profiles,
 *   non-disease-term profiles) using the prioritized gene list and publish
 *   outputs into the Druggability folder.
 *
 * Inputs:
 *   - prioritised_ch : tuple(runroot_path, Prioritised_genes.csv)
 *
 * Outputs:
 *   - results channel forwarded from the NON_DISEASE_TERM_DRUG_PROFILE process
 */
nextflow.enable.dsl=2

include { DRUG_TABLE }                    from '../modules/druggability/drug_table.nf'
include { GENE_DRUG_TARGET_PROFILES }     from '../modules/druggability/gene_drug_target_profiles.nf'
include { NON_DISEASE_TERM_DRUG_PROFILE } from '../modules/druggability/non_disease_term_drug_profile.nf'

/*
 * Input:
 *   prioritised_ch – tuples (runroot_path, path Prioritised_genes.csv)
 *
 * Emits:
 *   results – whatever NON_DISEASE_TERM_DRUG_PROFILE emits (same tuple)
 */

workflow DRUGGABILITY_FLOW {

  take:
    prioritised_ch

  main:
    def drug_in = prioritised_ch.map { runroot_path, prior_csv ->
      tuple( file(runroot_path), file(prior_csv) )
    }

    def drug_table_out      = DRUG_TABLE( drug_in )
    def gene_profiles_out   = GENE_DRUG_TARGET_PROFILES( drug_table_out )
    def nondisease_profiles = NON_DISEASE_TERM_DRUG_PROFILE( gene_profiles_out )

  emit:
    results = nondisease_profiles
}
