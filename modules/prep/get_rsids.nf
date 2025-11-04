/*
 * modules/prep/get_rsids.nf
 *
 * Purpose:
 *   Resolve RSIDs for per-trait lifted GWAS files by invoking the helper
 *   script that maps variant coordinates to rsids using a reference BIM/RDS.
 *
 * Inputs:
 *   - tuple(trait, runroot, lifted, script_get_rsids)
 *
 * Outputs:
 *   - <trait>_38_37_rsids.txt (symlinked into runroot)
 */
nextflow.enable.dsl=2

process GET_RSIDS {
    label 'light'
    label 'prep_batch5'
    tag { trait }

    input:
    tuple val(trait), path(runroot), path(lifted), path(script_get_rsids)

    output:
    tuple val(trait), path(runroot), path("runroot/${trait}_38_37_rsids.txt")

  script:
  """
    set -euo pipefail
    SCRIPT_PATH="\$(readlink -f "$script_get_rsids" 2>/dev/null || realpath "$script_get_rsids")"
    cd ${runroot}

    "\$SCRIPT_PATH" \
        --in "${lifted}" \
        --out "${params.output_path}/GWAS_Preprocessing/${trait}_38_37_rsids.txt" \
        --bim_rds "${params.databases}/1000genomes/EUR_phase3_autosomes.bim.rds"

    ln -sf "${params.output_path}/GWAS_Preprocessing/${trait}_38_37_rsids.txt" "./${trait}_38_37_rsids.txt"
  """

  stub:
  """
    mkdir -p runroot
    : > runroot/${trait}_38_37_rsids.txt
  """
}
