# HPC (cluster) run notes for nf-postgwas

Contents
--------
- Quick checklist to adapt `conf/hpc.config` to your site
- Example job scripts (SGE serial submit and MPI/ignite parallel)
- Apptainer (container) usage notes and best-practices
- Where to put site-specific overrides

Quick adaptation checklist
--------------------------
1. Pick the executor your site supports:
   - SGE: `process.executor = 'sge'`
   - SLURM: `process.executor = 'slurm'`
   - PBS/Torque: `process.executor = 'pbs'`
   - Ignite (MPI): omit the executor or set `process.executor = 'ignite'`
2. Update the queue/partition name in `conf/hpc.config` (`queue` for SGE,
   `partition` for SLURM).
3. Tune `clusterOptions` so Nextflow submits jobs with the right flags for
   your scheduler (examples are in `conf/hpc.config`).
4. Consider making a site-specific config (e.g., `conf/site.example.config`)
   and add it to `.gitignore` so users can copy/modify it safely.
5. Test with a tiny resource request (1 core, small memory, short walltime)
   before running full analyses.

Example job script — SGE (serial master job)
-------------------------------------------
Create a small script (e.g. `examples/run_nextflow_sge.sh`) and adapt the
resource flags to your queue/account.

```bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=24:0:0
#$ -l h_vmem=12G

# Option A: run Nextflow inside an Apptainer container (recommended)
# You may want to pull the container once on a compute node:
# apptainer pull postgwas-pipeline.sif docker://hlnicholls/postgwas-pipeline:latest

apptainer exec postgwas-pipeline.sif \
  nextflow -DXmx=1G -C nextflow.config run . -profile hpc -params-file params.yaml

# Option B: use the system-provided Nextflow module instead
# module load nextflow
# nextflow -DXmx=1G -C nextflow.config run . -profile hpc -params-file params.yaml
```

Example job script — MPI / Ignite (parallel)
-------------------------------------------
When you need many cores across nodes and your site supports Ignite+MPI,
request multiple nodes and run the master via `mpirun`.

```bash
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe parallel 96
#$ -l h_rt=240:0:0

module load openmpi
apptainer exec postgwas-pipeline.sif \
  mpirun --pernode nextflow run . -profile hpc -params-file params.yaml
```

Apptainer / container notes
---------------------------
- Pull large images on a compute node / interactive session to avoid long
  pulls inside short batch jobs.
- Example pull (once):
  `apptainer pull postgwas-pipeline.sif docker://hlnicholls/postgwas-pipeline:latest`
- If your site provides Apptainer-wrapped modules for common containers, you
  can `module load postgwas-pipeline` and use the module name instead of the
  `apptainer exec` command.
- Prefer `apptainer exec` when you want to run specific commands inside the
  container; `apptainer run` will execute the container's runscript.

Site-specific overrides
-----------------------
- Create a `conf/site.config` or `conf/<site>.config` (not committed) with
  your site-specific queue names, account flags, and any `clusterOptions`
  needed. Example:

```nextflow
profiles { site { includeConfig 'conf/hpc.config'; includeConfig 'conf/site.apocrita.config' } }
```

and a `conf/site.apocrita.config` (local) can contain only the lines you want
to change (queue name, account flags, etc.).
