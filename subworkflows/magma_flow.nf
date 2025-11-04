/*
 * subworkflows/magma_flow.nf
 *
 * Purpose:
 *   Run MAGMA gene-based analysis for each trait. This workflow
 *   prepares RSID inputs and invokes MAGMA, optionally skipping if
 *   results already exist.
 *
 * Inputs:
 *   - prep_out_ch : channel from PREP_FLOW with (trait, runroot, rsids_file)
 *
 * Outputs:
 *   - magma_out : channel with per-trait MAGMA outputs
 */
nextflow.enable.dsl=2

include { PREP_FOR_MAGMA } from '../modules/magma/prep_for_magma'
include { RUN_MAGMA }      from '../modules/magma/run_magma'

workflow MAGMA_FLOW {
  take:
    prep_out_ch   // tuple val(trait), val(runroot_path), path(rsids_file) from PREP_FLOW

  main:
    // Check if MAGMA results already exist for all traits  
    def allTraits = params.traits ?: []
  def magmaResultsDir = file("${params.output_path}/MAGMA/Results")
    def allResultsExist = false
    
    if (magmaResultsDir.exists() && allTraits.size() > 0) {
      // Check for both .genes.out and .genes.raw files for each trait
      allResultsExist = allTraits.every { trait ->
        def genesOut = file("${magmaResultsDir}/magma_${trait}.genes.out")
        def genesRaw = file("${magmaResultsDir}/magma_${trait}.genes.raw")
        genesOut.exists() && genesRaw.exists()
      }
      
      if (allResultsExist) {
        log.info "⊘ Skipping MAGMA analysis - results already exist for all ${allTraits.size()} traits"
      }
    }
    
    if (allResultsExist) {
      // Create empty channels - skip both PREP_FOR_MAGMA and RUN_MAGMA
      magma_out_ch = Channel.empty()
    } else {
      log.info "✓ Running MAGMA gene-based analysis (${allTraits.size()} traits)"
  // Pass the correct tuple from PREP_FLOW output to PREP_FOR_MAGMA
  all_traits_ch = prep_out_ch.map { trait, runroot, rsids_file -> tuple(runroot, trait, rsids_file, false) }

  // Run prep and MAGMA for each trait in parallel
  prepped_ch = PREP_FOR_MAGMA(all_traits_ch)
  magma_out_ch = RUN_MAGMA(prepped_ch)
    }

  emit:
    magma_out = magma_out_ch
}
