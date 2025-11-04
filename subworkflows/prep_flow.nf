/*
 * subworkflows/prep_flow.nf
 *
 * Purpose:
 *   Prepare per-trait runroot configurations and optionally run the
 *   PREP sequence (merge chr outputs, QC, liftover, get rsids). Emits
 *   a channel of (trait, runroot, rsids_txt) used by downstream flows.
 *
 * Inputs:
 *   - ch_tuples : list of (trait, sample_size, genome_build, merged_file|null)
 *
 * Outputs:
 *   - rsids channel : tuples (trait, runroot, rsids_txt)
 */
nextflow.enable.dsl = 2

include { LINK_PREMADE_GWAS }   from '../modules/prep/link_premade_gwas.nf'
include { CREATE_CONFIG_SHIMS } from '../modules/utils/make_config_shims.nf'
include { MERGE_CHRS }          from '../modules/prep/merge_chrs.nf'
include { INITIAL_QC }          from '../modules/prep/initial_qc.nf'
include { LIFTOVER_38_TO_37 }   from '../modules/prep/liftover_38to37.nf'
include { GET_RSIDS }           from '../modules/prep/get_rsids.nf'

workflow PREP_FLOW {

    take:
        ch_tuples
        skip_processing  // Boolean: if true, skip MERGE/QC/LIFTOVER/RSIDS but still create configs

    main:
        // 1) Always produce config shims; downstream relies on (trait, runroot)
        confs = CREATE_CONFIG_SHIMS(
            ch_tuples.map { trait, ss, genome_build, _ -> tuple(trait, ss, genome_build) }
        )
        // confs shape: (trait, runroot)

        if( !skip_processing ) {
            // 2) Split input by presence of premade merged file
            merged_ch_filt   = ch_tuples.filter { trait, _ss, _gb, merged_file -> merged_file != null }
            to_merge_ch_filt = ch_tuples.filter { trait, _ss, _gb, _merged_file -> _merged_file == null }

            // 3) Prepare joinable pairs keyed by trait
            confs_pairs      = confs.map             { trait, runroot -> tuple(trait, runroot) }
            premade_pairs    = merged_ch_filt.map    { trait, _ss, _gb, mf -> tuple(trait, file(mf)) }
            to_merge_pairs   = to_merge_ch_filt.map  { trait, _ss, _gb, _ -> tuple(trait, true) }

            // 4) Join each branch with its runroot from confs
            merged_linked = premade_pairs
                .join(confs_pairs)                       // (trait, mf, runroot)
                .map { trait, mf, runroot ->
                    tuple(trait, runroot, mf)
                } \
                | LINK_PREMADE_GWAS
            // emits: (trait, path("runroot/${trait}_regenie_allchr.txt"))

            merged_merged = to_merge_pairs
                .join(confs_pairs)                       // (trait, _flag, runroot)
                .map { trait, _flag, runroot ->
                    tuple(trait, runroot, file("${workflow.projectDir}/bin/prep/step1_merge_chrs.sh"))
                } \
                | MERGE_CHRS
            // emits: (trait, path("runroot/${trait}_regenie_allchr.txt"))

            // 5) Unify both sources of merged GWAS
            merged = merged_linked.mix(merged_merged)

            // 6) QC
            qc = merged.map { trait, merged_path ->
                // Use the runroot directory that contains merged_path
                // nextflow will stage runroot as a directory for INITIAL_QC
                tuple(trait, merged_path.parent, file("${workflow.projectDir}/bin/prep/step2_initial_QC.py"))
            } | INITIAL_QC
            // -> (trait, runroot, qc_merged_txt)

            // 7) Liftover
            lifted = qc.map { trait, runroot, qc_merged ->
                tuple(trait, runroot, qc_merged, file("${workflow.projectDir}/bin/prep/step3_LiftOver_38to37.py"))
            } | LIFTOVER_38_TO_37
            // -> (trait, runroot, lifted_txt)

            // 8) RSIDs
            rsids = lifted.map { trait, runroot, lifted_file ->
                tuple(trait, runroot, lifted_file, file("${workflow.projectDir}/bin/prep/step4_get_rsids.R"))
            } | GET_RSIDS
            // -> (trait, runroot, rsids_txt)
        }
        else {
            // Skip processing: surface existing files with the same output shape
            rsids = confs.map { trait, runroot ->
                def rsidsFile = file("${params.output_path}/GWAS_Preprocessing/${trait}_38_37_rsids.txt")
                tuple(trait, runroot, rsidsFile)
            }
        }

    emit:
        rsids
}
