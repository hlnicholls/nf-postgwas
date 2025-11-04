/*
 * modules/ld/collate_ld.nf
 *
 * Purpose:
 *   Collate per-locus LD files into a global set and invoke the Python
 *   collation script to produce a combined loci table with LD annotations.
 *
 * Inputs:
 *   - gate, tuple(runroot_path, loci_flag, script_collate_ld)
 *
 * Outputs:
 *   - all_traits_loci_38_with_ld.txt
 */
nextflow.enable.dsl = 2

process COLLATE_LD {
  label 'medium'
  tag "all_traits"

  input:
    val gate
    tuple val(runroot_path), path(loci_flag), path(script_collate_ld)

  output:
    path "all_traits_loci_38_with_ld.txt"

  script:
  """
  set -euo pipefail
  WORK_DIR="\$PWD"

  # Resolve script path before changing directories
  SCRIPT_PATH="\$(readlink -f "$script_collate_ld" 2>/dev/null || realpath "$script_collate_ld")"

  cd ${runroot_path}

  # Ensure the global LD directory exists and is populated (the Python script expects it)
  GLOBAL_LD_DIR="${params.output_path}/Loci_Preprocessing/Single_trait_LD"
  LOCAL_LD_DIR="Loci_Preprocessing/Single_trait_LD"
  mkdir -p "\${GLOBAL_LD_DIR}"

  # If global dir is empty but local has files, sync them across
  shopt -s nullglob
  local_ld=(\${LOCAL_LD_DIR}/*.ld)
  global_ld=(\${GLOBAL_LD_DIR}/*.ld)
  if [ "\${#global_ld[@]}" -eq 0 ] && [ "\${#local_ld[@]}" -gt 0 ]; then
    rsync -a --exclude='*.log' "\${LOCAL_LD_DIR}/" "\${GLOBAL_LD_DIR}/"
    # refresh listing
    global_ld=(\${GLOBAL_LD_DIR}/*.ld)
  fi

  if [ "\${#global_ld[@]}" -eq 0 ]; then
    echo "ERROR: No .ld files found for collation in \${GLOBAL_LD_DIR}" >&2
    exit 3
  fi

  # WARNING: May be computationally intensive on large LD sets.
  "\$SCRIPT_PATH"

  cp "${params.output_path}/Loci_Preprocessing/all_traits_loci_38_with_ld.txt" \
     "\$WORK_DIR/all_traits_loci_38_with_ld.txt"
  """
}

workflow COLLATE_LD_FLOW {
  take:
    ld_done_gate    // a single true once all CALC_LD tasks finish
    runroot_path    // representative runroot (Channel emitting a single path)
    loci_flag       // path('loci_done.flag') singleton

  main:
    rf = runroot_path.combine(loci_flag).map { rr, flag ->
      tuple(rr, flag, file("${workflow.projectDir}/bin/ld/collate_LD.py"))
    }
    out_file = COLLATE_LD(ld_done_gate, rf)

  emit:
    out_file
}
