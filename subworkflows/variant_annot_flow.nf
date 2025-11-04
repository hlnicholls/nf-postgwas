/*
 * subworkflows/variant_annot_flow.nf
 *
 * Purpose:
 *   Run variant annotation submodules (CADD, pCHiC, RegulomeDB, VCF, VEP)
 *   and combine outputs into a single annotated variants table.
 *
 * Inputs:
 *   - ld_out_ch : LD-annotated variant file channel
 *
 * Outputs:
 *   - variant_annotations : combined variant annotation channel
 */
nextflow.enable.dsl=2

include { VARIANT_ANNOT_CADD }        from '../modules/variant_annot/cadd_from_source'
include { VARIANT_ANNOT_PCHIC }       from '../modules/variant_annot/pchic'
include { VARIANT_ANNOT_REGULOME }    from '../modules/variant_annot/regulome'
include { VARIANT_ANNOT_TO_VCF }      from '../modules/variant_annot/to_vcf'
include { VARIANT_ANNOT_COMBINE }     from '../modules/variant_annot/combine_annotations'

workflow VARIANT_ANNOT_FLOW {
  take:
    // From LD_FLOW: tuple val(runroot_path), path("all_traits_loci_38_with_ld_genes.txt")
    ld_out_ch

  main:
  // Check if final combined annotation file already exists
  def annotDir = file("${params.output_path}/Variant_annotation")
  def finalOutputFile = file("${annotDir}/all_variant_annotations.csv")
    def finalOutputExists = finalOutputFile.exists() && finalOutputFile.size() > 0

    if (finalOutputExists) {
      log.info "⊘ Skipping variant annotation - final combined results already exist"
      // Emit the existing file so downstream processes can still run
      combined_annot_ch = Channel.fromPath(finalOutputFile)
    } else {
        // VEP modules are not present in `main` (developed on `new-feats`).
        // Run the remaining annotators and combine results; VEP is skipped.
        log.info "✓ Running variant annotation workflow (CADD → pCHiC → RegulomeDB → VCF → Combine; VEP skipped)"

      // Step 1–4: Individual annotators (each must emit: tuple val(runroot_path), path(<file>))
      cadd_out_ch     = VARIANT_ANNOT_CADD(ld_out_ch)
      pchic_out_ch    = VARIANT_ANNOT_PCHIC(ld_out_ch)
      regulome_out_ch = VARIANT_ANNOT_REGULOME(ld_out_ch)
      vcf_out_ch      = VARIANT_ANNOT_TO_VCF(ld_out_ch)

      // Synchronise non-VEP paths; pass ONLY runroot to the combiner
      ready_for_combine = cadd_out_ch
        .join(pchic_out_ch)
        .join(regulome_out_ch)
        .join(vcf_out_ch)
        .map { runroot, _cadd, _pchic, _reg, _vcf -> runroot }

      // Step 7: Combine (expects: val runroot_path), writes all_variant_annotations.csv
      combined_annot_ch = VARIANT_ANNOT_COMBINE(ready_for_combine)
    }

  emit:
    variant_annotations = combined_annot_ch
}
