/*
 * subworkflows/fine_mapping_flow.nf
 *
 * Purpose:
 *   Run the configured fine-mapping method (default: Wakefield), annotate
 *   credible sets and optionally annotate with GTEx colocalisation results.
 *
 * Inputs:
 *   - ld_out_ch
 *   - variant_annot_out_ch
 *   - gtex_coloc_out_ch
 *
 * Outputs:
 *   - fine_mapping_results
 *   - gwas_catalog_results
 */
nextflow.enable.dsl=2

include { WAKEFIELD_CREDIBLE_SETS }        from '../modules/fine_mapping/wakefield'
include { ANNOTATE_CREDSETS_WAKEFIELD }    from '../modules/fine_mapping/annotate_credsets'
include { ANNOTATE_CREDSETS_WITH_GTEX }    from '../modules/fine_mapping/annotate_credsets'
include { GWAS_CATALOG_QUERY }             from '../modules/fine_mapping/gwas_catalog_query'

workflow FINE_MAPPING_FLOW {
  take:
    ld_out_ch             // from LD_FLOW - tuple val(runroot), path(ld_genes)  
    variant_annot_out_ch  // from VARIANT_ANNOT_FLOW - path(annotations) - dependency signal
    gtex_coloc_out_ch     // from COLOC_FLOW - for GTEx annotation step

  main:
    // Determine which fine-mapping method to run based on params
    def fineMappingMethod = params.fine_mapping_method ?: 'wakefield'
    
    log.info "✓ Running fine-mapping analysis using method: ${fineMappingMethod}"
    
    // Ensure variant annotation completes before fine-mapping
    // Combine ld_out_ch with variant_annot completion signal
    ready_ch = ld_out_ch.combine(variant_annot_out_ch.collect()).map { runroot, ld_genes, _annot -> 
      tuple(runroot, ld_genes)
    }
    
    // Run the selected fine-mapping method
    if (fineMappingMethod.toLowerCase() == 'wakefield') {
        wakefield_ch = WAKEFIELD_CREDIBLE_SETS(ready_ch)
        
        // Step 2: Annotate credible sets with variant annotations (method-specific)
        annotated_ch = ANNOTATE_CREDSETS_WAKEFIELD(wakefield_ch)
        
        // Step 3: Annotate with GTEx coloc results (method-agnostic)
        // The ANNOTATE_CREDSETS_WITH_GTEX process will handle missing files gracefully
        // by only annotating traits that have GTEx results available
        log.info "  → Will check for GTEx coloc results and annotate if available"
        gtex_marker = gtex_coloc_out_ch.collect()
        fine_mapping_ch = ANNOTATE_CREDSETS_WITH_GTEX(annotated_ch, gtex_marker)
      } 
      // Future methods can be added here:
      // else if (fineMappingMethod.toLowerCase() == 'susie') {
      //   susie_ch = SUSIE_FINE_MAPPING(ready_ch)
      //   annotated_ch = ANNOTATE_CREDSETS_SUSIE(susie_ch)  // Different annotation script
      //   gtex_marker = gtex_coloc_out_ch.collect()
      //   fine_mapping_ch = ANNOTATE_CREDSETS_WITH_GTEX(annotated_ch, gtex_marker)
      // }
      // else if (fineMappingMethod.toLowerCase() == 'finemap') {
      //   finemap_ch = FINEMAP_ANALYSIS(ready_ch)
      //   annotated_ch = ANNOTATE_CREDSETS_FINEMAP(finemap_ch)  // Different annotation script
      //   gtex_marker = gtex_coloc_out_ch.collect()
      //   fine_mapping_ch = ANNOTATE_CREDSETS_WITH_GTEX(annotated_ch, gtex_marker)
      // }
      else {
        error "Unknown fine-mapping method: ${fineMappingMethod}. Supported methods: wakefield"
      }
    
    // Step 4: Query GWAS Catalog for pleiotropy analysis (method-agnostic, runs after annotation)
    // This step uses the standardized annotated credible sets output regardless of fine-mapping method
    log.info "  → Querying GWAS Catalog for credible set variants"
    gwas_catalog_ch = GWAS_CATALOG_QUERY(annotated_ch)

  emit:
    fine_mapping_results = fine_mapping_ch
    gwas_catalog_results = gwas_catalog_ch
}
