nextflow.enable.dsl = 2

// This test pipeline runs only the CREATE_CONFIG_SHIMS process from the
// real pipeline. It creates a small runroot with config_R.R / config_python.py
// so downstream steps can be validated quickly.

include { CREATE_CONFIG_SHIMS } from '../../modules/utils/make_config_shims.nf'

// Copy produced config into a test-accessible reports folder.
process COPY_CONFIG {
  tag { trait }
  input:
  tuple val(trait), path(runroot_dir)

  script:
  """
  set -euo pipefail
  mkdir -p "${params.outdir}"
  cp "${runroot_dir}/config_R.R" "${params.outdir}/config_R_${trait}.R"
  """
}

workflow {
  // Provide a single trait tuple: (trait, sample_size, genome_build)
  def inCh = channel.value( tuple('smoketrait', 123, 'GRCh38') )

  // Run the process that writes the runroot/config_*.{R,py,sh}
  def out = CREATE_CONFIG_SHIMS(inCh)

  // Run the copy step to place a file under params.outdir for test assertions
  COPY_CONFIG(out)
}
