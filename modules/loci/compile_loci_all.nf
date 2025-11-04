/*
 * modules/loci/compile_loci_all.nf
 *
 * Purpose:
 *   Create the combined loci table across traits by invoking the R
 *   compile_loci.R script. Writes a 'loci_done.flag' sentinel when done.
 *
 * Inputs:
 *   - tuple pairs (list of [trait, rsids_path]) and script path
 *
 * Outputs:
 *   - loci_done.flag
 */
nextflow.enable.dsl = 2

process COMPILE_LOCI_ALL {
  tag "all_traits"
  label 'light'

  input:
  tuple val(pairs), path(script_compile_loci)   // List of [trait, rsids_path]

  output:
  path "loci_done.flag"

  script:
  def traitNames = pairs.collect{ it[0] as String }
  def sampleSizes = (params.sample_sizes ?: []).collect{ it as Integer }
  def ssVec = (0..<traitNames.size()).collect{ i -> (sampleSizes && i < sampleSizes.size()) ? sampleSizes[i] : 0 }
  def traitsR = traitNames.collect{ "'${it.replace("'", "\\'")}'" }.join(', ')
  def ssR     = ssVec.join(', ')

  """
  set -euo pipefail

  # Ensure output directory exists for loci compilation results
  mkdir -p "${params.output_path}/Loci_Preprocessing"

  # Create config_R.R for compile_loci.R
  cat > config_R.R << 'EOF'
traits <- c(${traitsR})
sample_sizes <- c(${ssR})
home_path <- getwd()
input_path <- file.path(getwd(), "GWAS_Input")
output_path <- "${params.output_path}"
databases <- "${params.databases}"
# gwas_path: directory containing per-trait *_38_37_rsids.txt from PREP flow
gwas_path <- file.path(output_path, "GWAS_Preprocessing")
EOF

  # Create minimal dirs expected by the R code
  mkdir -p GWAS_Input

  # Run compile_loci.R using the staged script path
  Rscript "$script_compile_loci"

  echo "loci_done" > loci_done.flag
  """
}
