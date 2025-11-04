/*
 * subworkflows/pops_flow.nf
 *
 * Purpose:
 *   Orchestrate PoPS execution: run per-trait PoPS predictions, clean results,
 *   and run per-locus aggregation. Skips work if final outputs already exist.
 *
 * Inputs:
 *   - _magma_out_ch : MAGMA completion channel (may be unused)
 *   - close_pvals_ch: loci grouping/close p-values flag
 *   - ld_out_ch     : LD annotated output channel
 *
 * Outputs:
 *   - pops_results : per-locus PoPS aggregated results
 */
nextflow.enable.dsl=2

include { RUN_POPS }         from '../modules/pops/run_pops'
include { CLEAN_POPS_RESULTS } from '../modules/pops/clean_results'
include { POPS_PER_LOCI }    from '../modules/pops/pops_per_loci'

workflow POPS_FLOW {
  take:
    _magma_out_ch         // from MAGMA_FLOW - ensures MAGMA completes before POPS (not directly used)
    close_pvals_ch        // from LOCI_GROUPS_FLOW - path("All_loci_blocks_with_close_pvals.csv")
    ld_out_ch             // from LD_FLOW - tuple val(runroot_path), path("all_traits_loci_38_with_ld_genes.txt")

  main:
  // Check if PoPS results already exist
  def popsDir = file("${params.output_path}/Pops")
    def finalOutputsExist = false
    
    if (popsDir.exists()) {
      def expectedFiles = [
        "pops_top_genes_per_locus.txt",
        "gwas_all_loci_top_pops_genes.txt",
        "pops_results_all_features_cleaned.csv"
      ]
      def existingFiles = popsDir.listFiles()*.name
      finalOutputsExist = expectedFiles.every { existingFiles.contains(it) }
      
      if (finalOutputsExist) {
        log.info "⊘ Skipping PoPS analysis - final results already exist"
      }
    }
    
    if (finalOutputsExist) {
      // Create empty channels
      pops_preds_ch = Channel.empty()
      cleaned_ch = Channel.empty()
      pops_per_loci_ch = Channel.empty()
    } else {
      log.info "✓ Running PoPS (Polygenic Priority Scores) analysis"
      
      // Get runroot from LD output (which has runroot)
  runroot_val = ld_out_ch.map { runroot, _ld_file -> runroot }
      
      // Wait for MAGMA to complete if it's running, otherwise proceed immediately
      // Count MAGMA outputs - if none (MAGMA was skipped), create a dummy completion signal
      magma_completion = _magma_out_ch
        .toList()
        .map { list -> list.isEmpty() ? 'magma_skipped' : 'magma_done' }
      
      // Create channel for each trait, gated by MAGMA completion signal
      all_traits_ch = magma_completion
        .combine(Channel.fromList(params.traits ?: []))
        .combine(runroot_val)
        .map { _signal, trait, runroot -> tuple(runroot, trait, false) }
      
      // Run PoPS per trait in parallel
      pops_preds_ch = RUN_POPS(all_traits_ch)
      
      // Collect all .preds files and group by runroot for cleaning step
      grouped_preds = pops_preds_ch
        .map { runroot, _trait, preds -> tuple(runroot, preds) }
        .groupTuple()
      
      // Clean results (processes all traits together)
      cleaned_ch = CLEAN_POPS_RESULTS(grouped_preds)
      
      // Get loci blocks file (wait for LOCI_GROUPS_FLOW to complete first)
      // Combine close_pvals signal with cleaned results to get runroot
      loci_blocks_ch = close_pvals_ch
        .combine(cleaned_ch)
        .map { _close_pvals_file, runroot, _cleaned ->
          tuple(runroot, file("${params.output_path}/Loci_Preprocessing/All_loci_blocks.csv"))
        }
      
      ld_genes_ch = ld_out_ch  // Already has correct format
      
      // Run per-loci analysis
      pops_per_loci_ch = POPS_PER_LOCI(
        cleaned_ch,
        loci_blocks_ch,
        ld_genes_ch
      )
    }

  emit:
    pops_results = pops_per_loci_ch
}
