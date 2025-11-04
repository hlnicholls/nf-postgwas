/*
 * subworkflows/enrichment_flow.nf
 *
 * Purpose:
 *   Run EnrichR gene annotation when outputs are not already present. Emits
 *   a channel with EnrichR outputs for downstream prioritisation.
 *
 * Inputs:
 *   - ldOutCh : tuples (runroot, loci_with_ld_genes)
 *
 * Outputs:
 *   - enrichr_out : channel containing EnrichR annotation results
 */
nextflow.enable.dsl=2

include { ENRICHR_ANNOTATION } from '../modules/enrichment/enrichr_annotation'

workflow ENRICHMENT_FLOW {
  take:
    ldOutCh   // tuple val(runroot_path), path("all_traits_loci_38_with_ld_genes.txt") from LD_FLOW

  main:
  // Check if EnrichR output already exists
  def enrichrDir = file("${params.output_path}/Enrichment")
    def existingFiles = enrichrDir.exists() ? 
                        enrichrDir.listFiles().findAll { it.name.matches(/All_genes_annotated_with_EnrichR_.*\.csv/) } : 
                        []
    
    if (existingFiles) {
      log.info "⊘ Skipping EnrichR annotation - output already exists: ${existingFiles*.name.join(', ')}"
      // Create empty channel for output
      enrichr_out_ch = Channel.empty()
    } else {
      log.info "✓ Running EnrichR gene annotation"
      // Prepare input for ENRICHR_ANNOTATION
      // Combine LD output with the enrichr script
      enrichr_input = ldOutCh
        .map { runroot, loci_ld ->
          tuple(runroot, loci_ld,
                file("${workflow.projectDir}/bin/genes/enrichr_gene_annotation.R"))
        }
      
      enrichr_out_ch = ENRICHR_ANNOTATION(enrichr_input)
    }

  emit:
    enrichr_out = enrichr_out_ch
}
