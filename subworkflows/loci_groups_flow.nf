/*
 * subworkflows/loci_groups_flow.nf
 *
 * Purpose:
 *   Group nearby loci into blocks and identify lead SNPs with close p-values
 *   in other traits to support per-locus aggregation and downstream analyses.
 *
 * Inputs:
 *   - loci_with_ld_genes_ch : tuple(val(runroot_path), path(loci_with_ld_genes))
 *
 * Outputs:
 *   - loci_blocks
 *   - close_pvals
 */
nextflow.enable.dsl=2

include { LOCI_GROUPS }           from '../modules/loci_groups/loci_groups'
include { IDENTIFY_CLOSE_PVALS }  from '../modules/loci_groups/identify_close_pvals'

workflow LOCI_GROUPS_FLOW {
  take:
    loci_with_ld_genes_ch  // tuple val(runroot_path), path(loci_with_ld_genes)
  
  main:
    // Step 1: Create loci groups based on 1MB distance and LD
    loci_blocks_ch = LOCI_GROUPS(loci_with_ld_genes_ch)
    
    // Step 2: Identify leads with close p-values in other traits
    close_pvals_ch = IDENTIFY_CLOSE_PVALS(loci_blocks_ch)
  
  emit:
    loci_blocks = loci_blocks_ch
    close_pvals = close_pvals_ch
}
