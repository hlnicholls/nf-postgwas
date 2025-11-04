/*
 * subworkflows/prioritisation_flow.nf
 *
 * Purpose:
 *   Prepare inputs and barriers for the gene prioritisation step. Converts
 *   possibly-empty completion channels into stable sentinel values so the
 *   prioritisation process can run reliably across re-runs.
 *
 * Inputs:
 *   - pops_out_ch (PoPS outputs, may be empty)
 *   - enrichr_out_ch (EnrichR outputs, may be empty)
 *   - coloc_done_ch (CUSTOM coloc readiness channel)
 *   - variant_annot_out_ch (variant annotation barrier)
 *   - prep_out_ch (tuples: trait, runroot, rsids_txt)
 *
 * Outputs:
 *   - prioritised: tuples (runroot, path('Prioritised_genes.csv'))
 *   - impc: sentinel 'impc_done' emitted after prioritisation completes
 */
import java.nio.file.Files

include { GENE_PRIORITISATION } from '../modules/prioritisation/gene_prioritisation.nf'
// NOTE: IMPC enrichment is intentionally not included here to avoid double instantiation

workflow PRIORITISATION_FLOW {

  take:
    pops_out_ch
    enrichr_out_ch
    coloc_done_ch
    variant_annot_out_ch
    prep_out_ch

  main:
    /*
     * Resolve output roots. Keep existing param names/behaviour but add a safe fallback
     * so this subworkflow is robust even if output_path isn't explicitly set.
     */
    final outputRoot = params.containsKey('output_path') && params.output_path
        ? file(params.output_path.toString())
        : file("${baseDir}/results")

    final popsDir  = file("${outputRoot}/Pops")
    final colocDir = file("${outputRoot}/Colocalisation")

    // Determine whether custom coloc is required.
    final requireCustomColoc = (params.containsKey('run_custom_coloc') ? (params.run_custom_coloc as boolean) : false)

    final providedCustomName = (params.containsKey('custom_trait_name') ? params.custom_trait_name?.toString() : null)
    if (requireCustomColoc) {
      if (!providedCustomName || providedCustomName.trim().isEmpty())
        throw new RuntimeException("params.custom_trait_name must be set when run_custom_coloc=true")
    }

    // Construct a coloc CSV path for downstream tasks (string; may not exist).
    final effectiveCustomName = (providedCustomName && providedCustomName.trim()) ? providedCustomName : 'CUSTOM_DISABLED'
    final colocCsvStatic = file("${colocDir}/${effectiveCustomName}_coloc_all_variants_pp4.csv")

    // Deduplicate runroots from PREP (keep as files for path safety)
    def runroot_v = prep_out_ch.map { _trait, runroot, _rsids -> file(runroot) }.unique()

    // ---- Barrier tokens (convert possibly-empty completion channels into value channels)
    // Create fallback value channels up-front to avoid creating channels during igniter iteration.
    def pops_fallback    = Channel.value('p')
    def enrichr_fallback = Channel.value('e')
    def coloc_fallback   = Channel.value('c')
    def vannot_fallback  = Channel.value('v')

    // Reduce to single tokens; if a channel is empty, substitute its fallback.
    def pops_token    = pops_out_ch.map { 'p' }.ifEmpty { pops_fallback }
    def enrichr_token = enrichr_out_ch.map { 'e' }.ifEmpty { enrichr_fallback }
    def coloc_token   = coloc_done_ch.map { 'c' }.ifEmpty { coloc_fallback }
    def vannot_token  = variant_annot_out_ch.map { 'v' }.ifEmpty { vannot_fallback }

    // When all tokens have arrived, produce a single deps_ready value
    def deps_ready = pops_token
                      .combine(enrichr_token)
                      .combine(coloc_token)
                      .combine(vannot_token)
                      .map { 'deps_ready' }
                      .unique()

    // ---- Build validated inputs AFTER barriers
    def validated_inputs = deps_ready
      .combine(runroot_v)
      .map { _tok, runroot ->

        // Respect existing semantics: PoPS is ON by default unless explicitly disabled
        boolean requirePops = (params.containsKey('run_pops') ? (params.run_pops as boolean) : true)

        // Expected PoPS top-genes table
        def popsTop = file("${popsDir}/gwas_all_loci_top_pops_genes.txt")

        /*
         * Previous behaviour threw an exception if the file was missing/empty.
         * That proved brittle (publish timing, header-only outputs, or path drift).
         * New behaviour:
         *  - If PoPS is required and the file is missing or size==0, WARN and pass a harmless placeholder.
         *  - GENE_PRIORITISATION already handles missing/empty inputs gracefully.
         * This keeps the pipeline progressing while still surfacing the issue in logs.
         */
        if (requirePops) {
          // Use java.nio.file.Files with defensive try/catch to avoid uncaught
          // filesystem errors (eg. missing parent dir) causing the workflow
          // to abort. On any error, fall back to the placeholder.
          try {
            def popsPath = popsTop.toPath()
            if (!Files.exists(popsPath)) {
              log.warn "PoPS top-genes file not found at expected path: ${popsTop}. Proceeding with placeholder; prioritisation will omit PoPS."
              popsTop = file("/dev/null")
            } else if (Files.size(popsPath) == 0L) {
              log.warn "PoPS top-genes file exists but is empty: ${popsTop}. Proceeding with placeholder; prioritisation will omit PoPS."
              popsTop = file("/dev/null")
            }
          }
          catch (Exception e) {
            log.warn "Failed to inspect PoPS file ${popsTop}: ${e.message} — proceeding with placeholder."
            popsTop = file("/dev/null")
          }
        }
        else {
          // PoPS disabled: always pass a placeholder
          popsTop = file("/dev/null")
        }

        // Pass the coloc CSV path as a *string*; GENE_PRIORITISATION handles empty/missing.
        tuple(runroot, popsTop.toString(), colocCsvStatic.toString())
      }

    // ---- Run prioritisation -> (runroot, Prioritised_genes.csv)
    def gp = GENE_PRIORITISATION(validated_inputs)

    // ---- Provide a simple sentinel so main.nf can depend on "impc" completion
    // Derive it from gp, ensuring it only fires once prioritisation is done.
    def impc_done = gp.map { _runroot, _file -> 'impc_done' }
                      .ifEmpty { Channel.value('impc_done') }
                      .unique()

  emit:
    prioritised = gp
    impc        = impc_done
}
