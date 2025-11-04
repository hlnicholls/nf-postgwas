/*
 * modules/coloc/gtex_coloc.nf
 *
 * Purpose:
 *   Prepare GTEx eGene inputs (symlink if candidate files exist) and run
 *   GTEx colocalisation prep + coloc R scripts for supplied traits.
 *
 * Inputs:
 *   - tuple runroot_path, loci_with_ld_genes, script_prep_gtex, script_gtex_coloc
 *
 * Outputs:
 *   - per-trait *_eQTL_*.tsv / .txt files
 *   - gtex_coloc_completion.log
 */
nextflow.enable.dsl=2

process GTEX_COLOC {
  label 'heavy'
  tag "all_traits"

  publishDir "${params.output_path}/Colocalisation/eQTL",
             mode: 'copy',
             overwrite: true,
             pattern: "*_eQTL_*.{tsv,txt}"
  
  publishDir "${params.output_path}/Colocalisation",
             mode: 'copy',
             overwrite: true,
             pattern: "gtex_coloc_completion.log"

  input:
  tuple val(runroot_path), path(loci_with_ld_genes), path(script_prep_gtex), path(script_gtex_coloc)

  output:
  path "*_eQTL_*.{tsv,txt}"
  path "gtex_coloc_completion.log"

  script:
  // Check if GTEx coloc completion marker exists and is valid
  def colocDir = file("${params.output_path}/Colocalisation")
  def eqtlDir = file("${colocDir}/eQTL")
  def completionMarker = file("${colocDir}/gtex_coloc_completion.log")
  // Groovy-resolved GTEx DB path used in the shell script below. This is expanded by Groovy
  // at compile time so the shell sees an absolute path (e.g. /nf-postgwas/databases/gtex_v8_data)
  def dbGtex = "${params.databases}/gtex_v8_data"
  
  def allResultsExist = completionMarker.exists() && {
    def lines = completionMarker.readLines()
    def completedTraits = lines.find { it.startsWith("Completed traits:") }?.split(":")[1]?.trim()?.toInteger() ?: 0
    def totalTraits = lines.find { it.startsWith("Total traits:") }?.split(":")[1]?.trim()?.toInteger() ?: 0
    completedTraits == totalTraits && totalTraits == params.traits.size()
  }()
  
  if (allResultsExist) {
    """
    # GTEx coloc results already exist - copy to work directory
    cp "${eqtlDir}"/*_eQTL_*.tsv . 2>/dev/null || true
    cp "${eqtlDir}"/*_eQTL_*.txt . 2>/dev/null || true
    cp "${completionMarker}" gtex_coloc_completion.log
    """
  } else {
    """
    set -euo pipefail
    
    # Store the work directory
    WORK_DIR=\$PWD
    
    # Resolve script paths before changing directories
    SCRIPT_PREP_PATH="\$(readlink -f "$script_prep_gtex" 2>/dev/null || realpath "$script_prep_gtex")"
    SCRIPT_COLOC_PATH="\$(readlink -f "$script_gtex_coloc" 2>/dev/null || realpath "$script_gtex_coloc")"
    
    # Run the script from within the runroot
    cd ${runroot_path}
    
  # Ensure eQTL directory exists
  mkdir -p "${params.output_path}/Colocalisation/eQTL"
    
    # Some installations place files under gtex_eqtl_associations with names like
    # <TISSUE>.v8.EUR.allpairs.AllChr.txt.gz. Try to locate and symlink any candidate that
    # already contains a 'qval' column (i.e. gene-level results). If only allpairs exist we
    # do not attempt to convert them here — instead emit a clear advisory so the user can
    # provide proper eGene files or enable a separate conversion workflow.
    for tissue in ${params.gtex_tissues.join(' ')}; do
      EXPECTED="${dbGtex}/\$tissue.v8.egenes.txt.gz"
      if [ -f "\$EXPECTED" ]; then
        echo "Found expected GTEx egenes file: \$EXPECTED"
        continue
      fi

      # Look for candidate files under gtex_eqtl_associations (Groovy-expanded base path)
  CAND=`ls "${dbGtex}/gtex_eqtl_associations/\$tissue"* 2>/dev/null | head -n 1 || true`
      if [ -n "\${CAND}" ]; then
        echo "Found candidate GTEx file for \$tissue: \${CAND}"
        # Check for a header containing 'qval' (case-insensitive)
        if zcat "\${CAND}" 2>/dev/null | head -n1 | grep -qi 'qval'; then
          # Create a stable symlink at the expected path
          # Ensure the GTEx db directory exists (avoid nested command-substitution which can confuse the Groovy lexer)
          mkdir -p "${dbGtex}" >/dev/null 2>&1 || true
          ln -sf "\${CAND}" "${dbGtex}/\$tissue.v8.egenes.txt.gz"
          echo "Symlinked \${CAND} -> ${dbGtex}/\$tissue.v8.egenes.txt.gz"
        else
          echo "Candidate file \${CAND} does not contain a 'qval' column; cannot use as egenes file." >&2
          echo "Please provide GTEx gene-level eGenes (files named <TISSUE>.v8.egenes.txt.gz) in: ${dbGtex}" >&2
        fi
      else
        echo "No GTEx file found for tissue \$tissue (expected: ${dbGtex}/\$tissue.v8.egenes.txt.gz)" >&2
        echo "Please place GTEx eGenes under: ${dbGtex} or ${dbGtex}/gtex_eqtl_associations/ with an appropriate filename." >&2
      fi
    done

    # Step 1: Prep GTEx colocalisation (identify significant eQTLs)
    Rscript "\$SCRIPT_PREP_PATH"
    
    # Step 2: Run GTEx colocalisation analysis
    Rscript "\$SCRIPT_COLOC_PATH"
    
    # Create completion marker log
    echo "GTEx eQTL Colocalisation Completion Report" > "\$WORK_DIR/gtex_coloc_completion.log"
    echo "==========================================" >> "\$WORK_DIR/gtex_coloc_completion.log"
    echo "Timestamp: \$(date)" >> "\$WORK_DIR/gtex_coloc_completion.log"
    echo "Total traits: ${params.traits.size()}" >> "\$WORK_DIR/gtex_coloc_completion.log"
    
    # Count how many trait results were created
    COMPLETED=0
    for trait in ${params.traits.join(' ')}; do
      if [ -f "${params.output_path}/Colocalisation/eQTL/\${trait}_eQTL_COLOC.tsv" ]; then
        COMPLETED=\$((COMPLETED + 1))
        echo "  ✓ \${trait}: SUCCESS" >> "\$WORK_DIR/gtex_coloc_completion.log"
      else
        echo "  ✗ \${trait}: MISSING" >> "\$WORK_DIR/gtex_coloc_completion.log"
      fi
    done
    echo "Completed traits: \$COMPLETED" >> "\$WORK_DIR/gtex_coloc_completion.log"
    
  # Copy outputs back to work directory
  cp "${params.output_path}/Colocalisation/eQTL"/*_eQTL_*.tsv "\$WORK_DIR/" 2>/dev/null || true
  cp "${params.output_path}/Colocalisation/eQTL"/*_eQTL_*.txt "\$WORK_DIR/" 2>/dev/null || true
    """
  }
}
