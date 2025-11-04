/*
 * subworkflows/ld_flow.nf
 *
 * Purpose:
 *   Compute linkage-disequilibrium (LD) matrices for loci and produce
 *   annotated LD gene outputs for downstream modules.
 *
 * Inputs:
 *   - inTriples : tuples (trait, runroot, rsids_txt)
 *   - lociOut   : path flag indicating loci readiness
 *
 * Outputs:
 *   - ld_annotated : annotated LD results used by downstream flows
 */
nextflow.enable.dsl = 2

include { CALC_LD           } from '../modules/ld/calc_ld.nf'
include { COLLATE_LD_FLOW   } from '../modules/ld/collate_ld.nf'
include { ANNOTATE_LD_GENES } from '../modules/ld/annotate_genes.nf'

workflow LD_FLOW {
  take:
    inTriples  // (trait, runroot, rsids_txt)  -- across all traits
    lociOut    // path('loci_done.flag') singleton

  main:
    // Use a single representative tuple so CALC_LD runs ONCE globally
    rep = inTriples.take(1)

    // Gate by loci flag and conform to CALC_LD tuple signature
    ch_for_calc = rep.combine(lociOut)
      .map { trait, runroot, rsids_txt, flag ->
        tuple(trait as String, runroot, rsids_txt, flag, file("${workflow.projectDir}/bin/ld/calculate_LD.sh"))
      }

    // Single CALC_LD job that processes ALL loci (combined CSV)
    raw = CALC_LD(ch_for_calc)  // emits: (trait, runroot, runroot/<trait>_ld_raw.txt)

    // Barrier: single true when the lone CALC_LD finishes
    ld_done_gate = raw.map { true }

    // Reuse representative runroot (path) for downstream steps
    runroot_one = rep.map { it[1] }

    // Collate LD only after CALC_LD + loci flag
    col = COLLATE_LD_FLOW(ld_done_gate, runroot_one, lociOut)

    // Annotate LD genes - add script path
    anno_in = runroot_one.combine(col).map { rr, collated ->
      tuple(rr, collated, file("${workflow.projectDir}/bin/loci/annotate_genes.R"))
    }
    anno = ANNOTATE_LD_GENES(anno_in)

  emit:
    ld_annotated = anno
}
