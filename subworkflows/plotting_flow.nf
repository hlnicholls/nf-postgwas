/*
 * subworkflows/plotting_flow.nf
 *
 * Purpose:
 *   Run per-trait Manhattan/QC plotting and regional locus plots using LD
 *   annotations and compiled loci. Outputs plotting channels for downstream use.
 *
 * Inputs:
 *   - prepTriples : per-trait runroot metadata
 *   - ldAnnotated : LD annotation outputs
 *   - lociCompiled : compiled loci per trait
 *
 * Outputs:
 *   - manhattan : channel with Manhattan/QC plot artifacts
 *   - regional  : channel with regional locus plot artifacts
 */
nextflow.enable.dsl = 2

include { MANHATTAN_QC }   from '../modules/plots/manhattan_qc.nf'
include { REGIONAL_PLOTS } from '../modules/plots/regional_plots.nf'

workflow PLOT_FLOW {
  take:
    prepTriples
    ldAnnotated
    lociCompiled

  main:
    def qcPlotsR          = file("${projectDir}/bin/plots/QC_plots.R")
    def manhattanR        = file("${projectDir}/bin/plots/manhattan.R")
    def locuszoomR        = file("${projectDir}/bin/plots/locuszoom_plots.R")
    def makeRegionalPlotR = file("${projectDir}/bin/plots/makeRegionalPlot.R")

    // --- Manhattan + QC (per-trait; skippable inside the process)
    plotsIn = prepTriples.map { trait, runroot, rsids_file ->
      tuple(trait, runroot, rsids_file, qcPlotsR, manhattanR)
    }
    MANHATTAN_QC(plotsIn)

    // --- Regional plots (ALWAYS run; not skippable)
    // ldAnnotated emits (runroot, ld_file) — reduce to the LD file only
    ldFileOnly = ldAnnotated.map { rr, ld_file -> ld_file }

      regionalIn = prepTriples
        .combine(ldFileOnly)                           // (trait, runroot, rsids_file, ld_file)
        .map { trait, runroot, rsids_file, ld_file ->  // correct arity
          tuple(runroot, ld_file, locuszoomR, makeRegionalPlotR)   // (runroot, ld_file, script_locuszoom, script_makeRegionalPlot)
        }
      REGIONAL_PLOTS(regionalIn)

  emit:
    manhattan = MANHATTAN_QC.out
    regional  = REGIONAL_PLOTS.out
}
