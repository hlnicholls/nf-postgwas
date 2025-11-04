/*
 * modules/fine_mapping/annotate_credsets.nf
 *
 * Purpose:
 *   Annotate wakefield credible sets with variant-level annotations and
 *   optionally augment with GTEx coloc results when available.
 *
 * Inputs:
 *   - tuple(runroot_path, wakefield_results)
 *   - optional GTEx coloc marker (dependency)
 *
 * Outputs:
 *   - Annotated_CredibleSets_and_LD08.csv
 *   - Annotated_credsets_LD08_GTExColoc.csv (when GTEx available)
 */
nextflow.enable.dsl=2

process ANNOTATE_CREDSETS_WAKEFIELD {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Credible_set_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "Annotated_CredibleSets_and_LD08.csv"

  input:
  tuple val(runroot_path), path(wakefield_results)

  output:
  tuple val(runroot_path), path("Annotated_CredibleSets_and_LD08.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/fine_mapping/credible_set_annotation_wakefield.py"
  """
  set -euo pipefail

  WORK_DIR=\$PWD

  cd ${runroot_path}

  # Ensure output directory exists
  mkdir -p "${params.output_path}/Credible_set_annotation"

  # Use shared config (created by CREATE_CONFIG_SHIMS)
  python ${script_path}

  # Stage output back to work dir
  cp "${params.output_path}/Credible_set_annotation/Annotated_CredibleSets_and_LD08.csv" "\$WORK_DIR/"
  """
}

process ANNOTATE_CREDSETS_WITH_GTEX {
  label 'light'
  tag "all_traits"

  publishDir "${params.output_path}/Credible_set_annotation",
             mode: 'copy',
             overwrite: true,
             pattern: "Annotated_credsets_LD08_GTExColoc.csv"

  input:
  tuple val(runroot_path), path(annotated_credsets)
  path gtex_coloc_marker  // Just a dependency signal

  output:
  tuple val(runroot_path), path("Annotated_credsets_LD08_GTExColoc.csv")

  script:
  def script_path = "${workflow.projectDir}/bin/credsets/annotate_credset_with_GTEx_coloc.R"
  def outputFile = file("${params.output_path}/Credible_set_annotation/Annotated_credsets_LD08_GTExColoc.csv")
  
  if (outputFile.exists() && outputFile.size() > 0) {
    """
    echo "⊘ GTEx-annotated credible sets already exist, copying existing file..."
    cp "${outputFile}" Annotated_credsets_LD08_GTExColoc.csv
    """
  } else {
    """
    set -euo pipefail

    WORK_DIR=\$PWD

    # Touch input to ensure it exists
    test -s "${annotated_credsets}" > /dev/null

    cd ${runroot_path}

    # Ensure output directory exists
  mkdir -p "${params.output_path}/Credible_set_annotation"
    
    # Check if GTEx coloc results exist
  GTEX_DIR="${params.output_path}/Colocalisation/eQTL"
    if [ -d "\$GTEX_DIR" ] && ls "\$GTEX_DIR"/*_eQTL_COLOC.tsv 1> /dev/null 2>&1; then
      # Count available results
      TOTAL_TRAITS=\$(echo "${params.traits.join(' ')}" | wc -w)
      COMPLETED_TRAITS=0
      
      for trait in ${params.traits.join(' ')}; do
        if [ -f "\$GTEX_DIR/\${trait}_eQTL_COLOC.tsv" ]; then
          COMPLETED_TRAITS=\$((COMPLETED_TRAITS + 1))
        fi
      done
      
      echo "✓ Found GTEx coloc results for \$COMPLETED_TRAITS/\$TOTAL_TRAITS traits"
      echo "  Running GTEx annotation (traits without results will be skipped in R script)"
      
      # Use shared config (created by CREATE_CONFIG_SHIMS)
      Rscript ${script_path}
    else
      echo "⊘ No GTEx coloc results found - creating output without GTEx annotations"
      # Just copy the input file as output
  cp "${annotated_credsets}" "${params.output_path}/Credible_set_annotation/Annotated_credsets_LD08_GTExColoc.csv"
    fi

    # Stage output back to work dir
  cp "${params.output_path}/Credible_set_annotation/Annotated_credsets_LD08_GTExColoc.csv" "\$WORK_DIR/"
    """
  }
}
