/*
 * subworkflows/coloc_flow.nf
 *
 * Purpose:
 *   Coordinate CUSTOM and GTEx colocalisation modules and expose a
 *   readiness sentinel for downstream prioritisation.
 *
 * Inputs:
 *   - loci_csv_for_coloc_ch : tuples (trait, runroot_path, <ignored third>)
 *   - ld_annotated          : tuples (val runroot_path, path loci_with_ld_genes)
 *
 * Notes:
 *   - Individual modules (CUSTOM/GTEx) handle their own skip-if-exists
 *     behavior; this subworkflow only wires and signals readiness.
 */
nextflow.enable.dsl=2

include { GTEX_COLOC }   from '../modules/coloc/gtex_coloc'
include { CUSTOM_COLOC } from '../modules/coloc/custom_coloc'

workflow COLOC_FLOW {

  take:
  loci_csv_for_coloc_ch
  ld_annotated

  main:
  // Canonical all-loci CSV written by LOCI_FLOW
  def ALL_LOCI_CANON = file("${params.output_path}/Loci_Preprocessing/All_loci_ungrouped.csv")

  // Script paths (staged by Nextflow)
  def CUSTOM_COLOC_SCRIPT      = file("${baseDir}/bin/coloc/custom_coloc.R")
  def CUSTOM_COLOC_PROCESS     = file("${baseDir}/bin/coloc/process_custom_coloc_output.R")

  def GTEX_PREP_SCRIPT         = file("${baseDir}/bin/coloc/prep_GTEx_coloc.R")
  def GTEX_COLOC_SCRIPT_FILE   = file("${baseDir}/bin/coloc/GTEx_coloc.R")

  // Resolve toggles (with safe defaults)
  def _run_custom = params.containsKey('run_custom_coloc') ? (params.run_custom_coloc as boolean) : true
  def _run_gtex   = params.containsKey('run_gtex_coloc')   ? (params.run_gtex_coloc as boolean)   : true

  // Map to runroot (ignore any all_loci path from upstream to avoid wrong staging)
  def runroot_ch = loci_csv_for_coloc_ch.map { trait, runroot, _ignored -> runroot.toString() }

  // ---- CUSTOM (optional) — only run if loci CSV exists and is non-empty
  def custom_out_ch = nextflow.Channel.empty()
  if (_run_custom) {
    log.info "  • CUSTOM coloc: ENABLED"
    def custom_input = runroot_ch.map { runroot ->
      tuple(runroot, ALL_LOCI_CANON, CUSTOM_COLOC_SCRIPT, CUSTOM_COLOC_PROCESS)
    }
    // Only run CUSTOM_COLOC if the loci CSV exists and is non-empty
    def custom_input_filtered = custom_input.filter { runroot, loci_csv, script, process_script ->
      loci_csv.exists() && loci_csv.size() > 0
    }
    custom_out_ch = CUSTOM_COLOC(custom_input_filtered)
  } else {
    log.info "⊘ CUSTOM coloc: DISABLED"
  }

  // ---- GTEx (toggleable)
  if (_run_gtex) {
    log.info "✓ GTEx coloc: ENABLED"
    // Pair ld_annotated with runroot_ch. Using `combine` pairs the
    // emitted items by position which is robust when the channel values
    // may have different types (Path vs String) but represent the same
    // logical runroot for the workflow. This prevents an empty join
    // when the types don't match exactly.
    // The ld_annotated channel emits tuples (runroot_path, loci_ld_path).
    // When combined with runroot_ch the map closure receives three args
    // (runroot_from_ld, loci_ld_path, runroot_str) due to Groovy's tuple
    // destructuring. Accept three parameters to avoid invocation errors.
    def gtex_input = ld_annotated.combine(runroot_ch).map { runroot_from_ld, loci_ld_path, _runroot_str ->
      tuple(runroot_from_ld.toString(), loci_ld_path, GTEX_PREP_SCRIPT, GTEX_COLOC_SCRIPT_FILE)
    }
    GTEX_COLOC(gtex_input)
  } else {
    log.info "⊘ GTEx coloc: DISABLED (params.run_gtex_coloc=false)"
  }

  // Readiness sentinel for CUSTOM coloc:
  // - If CUSTOM is enabled, wait for one emission from its output channel
  // - If CUSTOM is disabled, emit a ready value so downstream doesn't stall
  def custom_ready = _run_custom
                      ? custom_out_ch.map { 'custom_coloc_ready' }
                      : channel.value('custom_coloc_ready')

  emit:
  custom = custom_ready
}
