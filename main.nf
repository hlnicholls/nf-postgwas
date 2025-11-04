/*
 * main.nf
 *
 * Purpose:
 *   Top-level workflow orchestration for the nf-postgwas pipeline. Wires
 *   subworkflows, toggles optional modules, and exposes channels used by
 *   downstream flows.
 *
 * Notes:
 *   - This file declares the OPT map of runtime toggles derived from params.
 */
nextflow.enable.dsl = 2

include { PREP_FLOW }            from './subworkflows/prep_flow.nf'
include { LOCI_FLOW }            from './subworkflows/loci_flow.nf'
include { LD_FLOW }              from './subworkflows/ld_flow.nf'
include { LOCI_GROUPS_FLOW }     from './subworkflows/loci_groups_flow.nf'
include { VARIANT_ANNOT_FLOW }   from './subworkflows/variant_annot_flow.nf'
include { FINE_MAPPING_FLOW }    from './subworkflows/fine_mapping_flow.nf'
include { PLOT_FLOW }            from './subworkflows/plotting_flow.nf'
include { LDSC_FLOW }            from './subworkflows/ldsc_flow.nf'
include { COLOC_FLOW }           from './subworkflows/coloc_flow.nf'
include { ENRICHMENT_FLOW }      from './subworkflows/enrichment_flow.nf'
include { MAGMA_FLOW }           from './subworkflows/magma_flow.nf'
include { POPS_FLOW }            from './subworkflows/pops_flow.nf'
include { DRUGGABILITY_FLOW }    from './subworkflows/druggability_flow.nf'
include { PRIORITISATION_FLOW }  from './subworkflows/prioritisation_flow.nf'

workflow {

  // -------- Required params & defaults --------
  if ( !params.traits || params.traits.size() == 0 )
    error "params.traits must be set (e.g. --traits '[\"LVEDV\",\"LVESV\"]')"

  if ( params.sample_sizes && params.traits.size() != params.sample_sizes.size() )
    log.warn "params.traits (${params.traits.size()}) and params.sample_sizes (${params.sample_sizes.size()}) differ."

  params.merged_regenie_files = (params.merged_regenie_files ?: [:])
  params.genome_build = (params.genome_build ?: 'GRCh38')

  // -------- Option toggles --------
  // Note: Only PREP_FLOW is mandatory. Everything else can be toggled.
  def OPT = [
    // Core graph-building steps (now all togglable)
    run_loci              : (params.containsKey('run_loci')              ? (params.run_loci as boolean)              : true),
    run_ld                : (params.containsKey('run_ld')                ? (params.run_ld as boolean)                : true),
    run_loci_groups       : (params.containsKey('run_loci_groups')       ? (params.run_loci_groups as boolean)       : true),
    run_plot              : (params.containsKey('run_plot')              ? (params.run_plot as boolean)              : true),

    // Optional analytics
    skip_prep             : (params.containsKey('skip_prep')             ? (params.skip_prep as boolean)             : false),
    run_variant_annotation: (params.containsKey('run_variant_annotation')? (params.run_variant_annotation as boolean): false),
    run_magma             : (params.containsKey('run_magma')             ? (params.run_magma as boolean)             : false),
    run_ldsc              : (params.containsKey('run_ldsc')              ? (params.run_ldsc as boolean)              : false),
    run_druggability      : (params.containsKey('run_druggability')      ? (params.run_druggability as boolean)      : false),
    run_pops              : (params.containsKey('run_pops')              ? (params.run_pops as boolean)              : true),
    run_prioritisation    : (params.containsKey('run_prioritisation')    ? (params.run_prioritisation as boolean)    : true),
    run_coloc             : (params.containsKey('run_coloc')             ? (params.run_coloc as boolean)             : true),
    run_fine_mapping      : (params.containsKey('run_fine_mapping')      ? (params.run_fine_mapping as boolean)      : true),
    run_enrichment        : (params.containsKey('run_enrichment')        ? (params.run_enrichment as boolean)        : true),
    run_gtex_coloc        : (params.containsKey('run_gtex_coloc')        ? (params.run_gtex_coloc as boolean)        : true)
  ] as Map<String,Boolean>

  log.info "Options: " + OPT.collect { k, v -> "${k}=${v}" }.sort().join('  ')

  // -------- Seed channel --------
  def traitTriplesCh = Channel.fromList(
    params.traits
      .withIndex()
      .collect { t, i ->
        def ss = (params.sample_sizes && i < params.sample_sizes.size()) ? params.sample_sizes[i] : 0
        def merged_file = params.merged_regenie_files.containsKey(t) ? params.merged_regenie_files[t] : null
        tuple(t as String, ss as Integer, params.genome_build as String, merged_file)
      }
      .findAll { tup -> tup[0] && tup[0] != 'nan' }
  )

  // -------- Auto-skip PREP if outputs exist --------
  boolean skipPrepEffective = OPT.skip_prep
  if ( !skipPrepEffective ) {
    def allPrepOutputsExist = params.traits.every { trait ->
      def prepDir = "${params.output_path}/GWAS_Preprocessing"
      file("${prepDir}/${trait}_regenie_allchr.txt").exists() &&
      file("${prepDir}/${trait}_38_37.txt").exists() &&
      file("${prepDir}/${trait}_38_37_rsids.txt").exists()
    }
    if ( allPrepOutputsExist ) {
      log.info ""
      log.info "✓ All PREP_FLOW processing outputs found in ${params.output_path}/GWAS_Preprocessing"
      log.info "  Files checked per trait:"
      log.info "    - *_regenie_allchr.txt (MERGE_CHRS + INITIAL_QC)"
      log.info "    - *_38_37.txt (LIFTOVER_38_TO_37)"
      log.info "    - *_38_37_rsids.txt (GET_RSIDS)"
      log.info "  Skipping PREP_FLOW processing steps (configs will still be created)."
      log.info "  To force re-run, add --skip_prep false or delete the output files."
      skipPrepEffective = true
    }
  }

  // -------- PREP (mandatory) --------
  def prep = PREP_FLOW(traitTriplesCh, skipPrepEffective)
  def prep_rsids_ch = prep.rsids

  // Aliases
  def prep_for_loci_ch           = prep_rsids_ch
  def prep_for_ld_ch             = prep_rsids_ch
  def prep_for_plot_ch           = prep_rsids_ch
  def prep_for_magma_ch          = prep_rsids_ch
  def prep_for_ldsc_ch           = prep_rsids_ch
  // Some feature-specific aliases were removed or archived; other aliases remain.
  def prep_for_prioritisation_ch = prep_rsids_ch

  // -------- LOCI (opt) --------
  def lociDoneCh = Channel.empty()
  def compiledLociCh = Channel.empty()
  if ( OPT.run_loci ) {
    log.info "✓ Running LOCI discovery & compilation"
    (lociDoneCh, compiledLociCh) = LOCI_FLOW(prep_for_loci_ch)
  } else {
    log.info "⊘ Skipping LOCI (run_loci = false)"
  }
  def loci_csv_for_plot_ch  = compiledLociCh
  def loci_csv_for_coloc_ch = compiledLociCh

  // -------- LD (opt) --------
  def ld_annotated_ch = Channel.empty()
  if ( OPT.run_ld ) {
    if ( OPT.run_loci && lociDoneCh ) {
      log.info "✓ Running LD annotation"
      ld_annotated_ch = LD_FLOW(prep_for_ld_ch, lociDoneCh)
    } else {
      log.warn "⚠ LD requested but LOCI is disabled or no lociDoneCh available — LD requires loci. Skipping LD."
    }
  } else {
    log.info "⊘ Skipping LD (run_ld = false)"
  }
  def ld_for_variant_annot_ch = ld_annotated_ch
  def ld_for_loci_groups_ch   = ld_annotated_ch
  def ld_for_plot_ch          = ld_annotated_ch
  def ld_for_coloc_ch         = ld_annotated_ch
  def ld_for_enrich_ch        = ld_annotated_ch
  def ld_for_finemap_ch       = ld_annotated_ch

  // -------- VARIANT ANNOT (opt) --------
  def variantAnnotOutCh = Channel.empty()
  if ( OPT.run_variant_annotation ) {
    if ( OPT.run_ld && ld_for_variant_annot_ch ) {
      log.info "✓ Running variant annotation (CADD, pCHiC, RegulomeDB)"
      variantAnnotOutCh = VARIANT_ANNOT_FLOW(ld_for_variant_annot_ch)
    } else {
      log.warn "⚠ Variant annotation requested but LD is unavailable — requires LD outputs. Skipping variant annotation."
      variantAnnotOutCh = Channel.empty()
    }
  } else {
    log.info "⊘ Skipping variant annotation (run_variant_annotation = false)"
    variantAnnotOutCh = ld_for_variant_annot_ch ?: Channel.empty()
  }
  def variantForFineMapCh = variantAnnotOutCh
  def variantForPriorCh   = variantAnnotOutCh

  // -------- LOCI GROUPS (opt) --------
  def closePvalsOut = Channel.empty()
  if ( OPT.run_loci_groups ) {
    if ( OPT.run_ld && ld_for_loci_groups_ch ) {
      log.info "✓ Running loci grouping / block construction"
      def _lociBlocksOut
      (_lociBlocksOut, closePvalsOut) = LOCI_GROUPS_FLOW(ld_for_loci_groups_ch)
    } else {
      log.warn "⚠ Loci groups requested but LD is unavailable — requires LD outputs. Skipping loci groups."
    }
  } else {
    log.info "⊘ Skipping loci groups (run_loci_groups = false)"
  }

  // -------- MAGMA (opt) --------
  def magmaOutCh = Channel.empty()
  if ( OPT.run_magma ) {
    log.info "✓ Running MAGMA gene-based analysis"
    magmaOutCh = MAGMA_FLOW(prep_for_magma_ch)
  } else {
    log.info "⊘ Skipping MAGMA analysis (run_magma = false)"
  }

  // -------- PoPS (opt) --------
  def popsOutCh = Channel.empty()
  if ( OPT.run_pops ) {
    if ( OPT.run_magma ) {
      if ( OPT.run_loci_groups && closePvalsOut && (OPT.run_ld && ld_for_loci_groups_ch) ) {
        log.info "✓ Running PoPS gene prioritisation"
        popsOutCh = POPS_FLOW(magmaOutCh, closePvalsOut, ld_for_loci_groups_ch)
      } else {
        log.warn "⚠ PoPS requested but loci groups/LD inputs are unavailable — running PoPS with MAGMA only may not be supported. Skipping PoPS."
      }
    } else {
      log.warn "⚠ PoPS requested but MAGMA is disabled — requires MAGMA results. Skipping PoPS."
    }
  } else {
    log.info "⊘ Skipping PoPS analysis (run_pops = false)"
  }

  // -------- Plotting (opt) --------
  if ( OPT.run_plot ) {
    if ( OPT.run_ld && ld_for_plot_ch && OPT.run_loci && loci_csv_for_plot_ch ) {
      log.info "✓ Running plotting"
      PLOT_FLOW(prep_for_plot_ch, ld_for_plot_ch, loci_csv_for_plot_ch)
    } else {
      log.warn "⚠ Plotting requested but LD/LOCI data are unavailable — requires both. Skipping plotting."
    }
  } else {
    log.info "⊘ Skipping plotting (run_plot = false)"
  }

  // -------- LDSC (opt) --------
  if ( OPT.run_ldsc ) {
    log.info "✓ Running LDSC (LD Score Regression)"
    LDSC_FLOW(prep_for_ldsc_ch)
  } else {
    log.info "⊘ Skipping LDSC analysis (run_ldsc = false)"
  }

  // -------- Coloc (opt) --------
  def customColocReady = Channel.empty()
  if ( OPT.run_coloc ) {
    if ( OPT.run_ld && ld_for_coloc_ch && OPT.run_loci && loci_csv_for_coloc_ch ) {
      log.info "✓ Running colocalisation analyses"
      if ( OPT.run_gtex_coloc ) log.info "  • GTEx coloc: ENABLED"
      else                       log.info "  • GTEx coloc: DISABLED"
      def colocFlow = COLOC_FLOW(loci_csv_for_coloc_ch, ld_for_coloc_ch)
      customColocReady = colocFlow.custom
    } else {
      log.warn "⚠ Coloc requested but LD and/or LOCI are unavailable — requires both. Skipping coloc."
      customColocReady = Channel.value('custom_coloc_skipped')
    }
  } else {
    log.info "⊘ Skipping colocalisation (run_coloc = false)"
    customColocReady = Channel.value('custom_coloc_disabled')
  }

  // -------- Fine mapping (opt) --------
  def fineMappingOutCh = Channel.empty()
  def gwasCatalogOutCh = Channel.empty()
  if ( OPT.run_fine_mapping ) {
    if ( OPT.run_ld && ld_for_finemap_ch ) {
      log.info "✓ Running fine mapping (method: ${params.fine_mapping_method})"
      (fineMappingOutCh, gwasCatalogOutCh) =
        FINE_MAPPING_FLOW(ld_for_finemap_ch, variantForFineMapCh, Channel.empty())
    } else {
      log.warn "⚠ Fine mapping requested but LD is unavailable — requires LD outputs. Skipping fine mapping."
    }
  } else {
    log.info "⊘ Skipping fine mapping (run_fine_mapping = false)"
  }

  // -------- Enrichment (opt) --------
  def enrichrOutCh = Channel.empty()
  if ( OPT.run_enrichment ) {
    if ( OPT.run_ld && ld_for_enrich_ch ) {
      log.info "✓ Running gene enrichment (gProfiler, IMPC)"
      enrichrOutCh = ENRICHMENT_FLOW(ld_for_enrich_ch)
    } else {
      log.warn "⚠ Enrichment requested but LD is unavailable — requires LD outputs. Skipping enrichment."
    }
  } else {
    log.info "⊘ Skipping enrichment (run_enrichment = false)"
  }

  // -------- Prioritisation (opt) --------
  def prioritisedCh  = Channel.empty()
  def impcSentinelCh = Channel.empty()
  if ( OPT.run_prioritisation ) {
    def missingDeps = [
      OPT.run_pops       != true ? "run_pops"       : null,
      OPT.run_enrichment != true ? "run_enrichment" : null,
      OPT.run_coloc      != true ? "run_coloc"      : null
    ].findAll { dep -> dep }
    if (missingDeps) {
      log.warn "⚠ Some dependencies for prioritisation are disabled: ${missingDeps.join(', ')}"
    }

    if ( (OPT.run_pops ? popsOutCh : true) &&
         (OPT.run_enrichment ? enrichrOutCh : true) ) {
      log.info "✓ Running gene prioritisation"
      def prFlow = PRIORITISATION_FLOW(
        popsOutCh,
        enrichrOutCh,
        customColocReady,
        variantForPriorCh,
        prep_for_prioritisation_ch
      )
      prioritisedCh  = prFlow.prioritised
      impcSentinelCh = prFlow.impc
    } else {
      log.warn "⚠ Prioritisation requested but required inputs are unavailable. Skipping prioritisation."
    }
  } else {
    log.info "⊘ Skipping prioritisation (run_prioritisation = false)"
  }

  // -------- Druggability (opt) --------
  def druggabilityOutCh = Channel.empty()
  if ( OPT.run_druggability ) {
    if ( OPT.run_prioritisation ) {
      if ( prioritisedCh ) {
        log.info "✓ Druggability will run after prioritisation completes"
        druggabilityOutCh = DRUGGABILITY_FLOW(prioritisedCh)
      } else {
        log.warn "⚠ Druggability requested but prioritisation outputs are unavailable. Skipping druggability."
      }
    } else {
      log.warn "⚠ Druggability requested but prioritisation is disabled — enable run_prioritisation to proceed."
    }
  } else {
    log.info "⊘ Skipping druggability (run_druggability = false)"
  }
}
