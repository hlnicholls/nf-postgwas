/*
 * subworkflows/loci_flow.nf
 *
 * Purpose:
 *   Compile locus definitions across all traits into a canonical loci table
 *   used by downstream LD and coloc modules.
 *
 * Inputs:
 *   - prepTriplesCh: tuples (trait, runroot, rsids_path)
 *
 * Outputs:
 *   - loci_done : path flag when compilation completes
 *   - compiled  : per-trait compiled loci CSV paths
 */
nextflow.enable.dsl = 2

include { COMPILE_LOCI_ALL } from '../modules/loci/compile_loci_all.nf'

workflow LOCI_FLOW {
  take:
    // From PREP: (trait, runroot, rsids_path)
    prepTriplesCh

  main:
    // Build [trait, rsids_path] pairs for the compile step
    lociPairs = prepTriplesCh.map { trait, runroot, rsids_path -> [trait, rsids_path] }
    allPairs  = lociPairs.toList()
    allPairsWithScript = allPairs.map { pairs -> tuple(pairs, file("${workflow.projectDir}/bin/loci/compile_loci.R")) }

    // Run a single compilation across all traits; emits a done flag
    lociDone = COMPILE_LOCI_ALL(allPairsWithScript)

    // Also emit per-trait expected compiled loci CSV path
    compiledPerTrait = prepTriplesCh.map { trait, runroot, rsids_path ->
  tuple(trait, runroot, file("${runroot}/All_loci_ungrouped.csv"))
    }

  emit:
    loci_done = lociDone            // path('loci_done.flag')
    compiled  = compiledPerTrait    // (trait, runroot, Singletrait_all_loci.csv)
}
