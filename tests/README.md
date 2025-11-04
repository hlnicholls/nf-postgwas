Smoke tests for nf-postgwas
==========================

This folder contains a minimal smoke test that runs a small part of the
pipeline's early setup: the CREATE_CONFIG_SHIMS process (from
`modules/utils/make_config_shims.nf`). That process writes per-trait
runroot configuration shims (`config_R.R`, `config_python.py`,
`config_shell.sh`) which downstream scripts rely on.

How to run
----------

- Using nf-test (preferred):

  From repository root run either:

  - `./nf-test test`  # run all discovered tests (recommended)
  - `./nf-test test tests/smoke/main.nf.test`  # run the smoke test explicitly

  nf-test will consult its configuration and execute the test descriptor
  `tests/smoke/main.nf.test` which runs `tests/smoke/main.nf` under the
  `test` profile.

- Direct Nextflow fallback:

  If you prefer to run Nextflow directly (no nf-test), use:

  ```bash
  nextflow run tests/smoke/main.nf -profile test -params-file tests/params.yaml -w .nf-test/work
  ```

Notes
-----
- The `-params-file tests/params.yaml` argument loads small,
  test-specific params so the main config doesn't require large external
  directories.
- The test writes a small work directory under `.nf-test/` which is created at
   runtime; this path is typically ignored by `.gitignore` and safe to remove
   between runs.
