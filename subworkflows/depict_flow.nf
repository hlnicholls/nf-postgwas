/*
 * subworkflows/depict_flow.nf
 *
 * Purpose:
 *   Orchestrate DEPICT pathway enrichment: create per-trait configs,
 *   prepare DEPICT inputs, run DEPICT, and collect formatted outputs.
 *
 * Inputs:
 *   - prep_out_ch : configuration/runroot channel from PREP_FLOW
 *
 * Outputs:
 *   - depict_results : prepared DEPICT output tables
 */
nextflow.enable.dsl=2

include { DEPICT_MAKE_CONFIG }    from '../modules/depict/make_config'
include { DEPICT_CREATE_INPUT }   from '../modules/depict/create_input'
include { RUN_DEPICT }            from '../modules/depict/run_depict'
include { DEPICT_PREPARE_OUTPUT } from '../modules/depict/prepare_output'

workflow DEPICT_FLOW {
  take:
    prep_out_ch  // tuple val(runroot_path), path(runroot_dir) - from PREP_FLOW configs

  main:
    log.info "✓ Running DEPICT pathway enrichment analysis"
    
    // Step 1: Create DEPICT configuration files per phenotype
    log.info "  → Creating DEPICT configuration files"
    config_ch = DEPICT_MAKE_CONFIG(prep_out_ch)
    
    // Step 2: Create DEPICT input files from GWAS results
    log.info "  → Creating DEPICT input files"
    input_ch = DEPICT_CREATE_INPUT(config_ch)
    
    // Step 3: Run DEPICT analysis
    log.info "  → Running DEPICT enrichment analysis (this may take a while)"
    depict_results_ch = RUN_DEPICT(input_ch)
    
    // Step 4: Prepare DEPICT output tables
    log.info "  → Preparing DEPICT output tables"
    depict_output_ch = DEPICT_PREPARE_OUTPUT(depict_results_ch)

  emit:
    depict_results = depict_output_ch
}
