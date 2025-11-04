/*
 * subworkflows/ldsc_flow.nf
 *
 * Purpose:
 *   Prepare GWAS inputs for LDSC, run heritability estimation and
 *   genetic correlation analyses, and clean the logs for publishing.
 *
 * Inputs:
 *   - inTriples : per-trait tuples (trait, runroot, rsids_txt)
 *
 * Outputs:
 *   - ldsc_results : cleaned LDSC outputs
 */
nextflow.enable.dsl = 2

include { CLEAN_FOR_LDSC    } from '../modules/ldsc/clean_for_ldsc.nf'
include { RUN_HERITABILITY  } from '../modules/ldsc/run_heritability.nf'
include { RUN_GENETIC_CORR  } from '../modules/ldsc/run_genetic_corr.nf'
include { CLEAN_LDSC_LOGS   } from '../modules/ldsc/clean_ldsc_logs.nf'

workflow LDSC_FLOW {
  take:
    inTriples  // (trait, runroot, rsids_txt)

  main:
    // Step 1: Clean GWAS data for LDSC format (per trait)
    // Add script path to the input tuple
    cleanInput = inTriples.map { trait, runroot, rsids_txt ->
      tuple(trait, runroot, rsids_txt, file("${workflow.projectDir}/bin/ldsc/clean_for_LDSC.py"))
    }
    cleaned = CLEAN_FOR_LDSC(cleanInput)
    
    // Step 2: Run heritability estimation (per trait)
    h2 = RUN_HERITABILITY(cleaned)

    // Step 3: Run genetic correlation for each trait against all others
    allH2 = h2.toList()
    gencorr = allH2.flatMap { items ->
      if (items.size() < 2) {
        // Not enough traits to compute genetic correlations
        return []
      }
      // items is a list of [trait, runroot, h2log]
      def allTraitNames = items.collect { it[0] }
      def allTraitFiles = allTraitNames.collect { t -> "${params.output_path}/GWAS_Preprocessing/${t}_GWAS_37_corr.txt" }

      items.collect { thisTrait ->
        def trait = thisTrait[0]
        def runroot = thisTrait[1]
        def trait_file = "${params.output_path}/GWAS_Preprocessing/${trait}_GWAS_37_corr.txt"
        def other_trait_files = allTraitFiles.findAll { it != trait_file }.join(',')
        tuple(trait, runroot, trait_file, other_trait_files)
      }
    } | RUN_GENETIC_CORR

    // Step 4: Clean and format the log files
    // Add script path to input
    cleanLogsInput = gencorr.map { trait, runroot, gencorr_log ->
      tuple(trait, runroot, gencorr_log, file("${workflow.projectDir}/bin/ldsc/clean_log_files.R"))
    }
    cleaned_logs = CLEAN_LDSC_LOGS(cleanLogsInput)

  emit:
    ldsc_results = cleaned_logs
}
