# PostGWAS Analysis Pipeline
[![Website](https://img.shields.io/website?down_color=red&down_message=down&up_message=up&url=https://hlnicholls.github.io/nf-postgwas/)](https://hlnicholls.github.io/nf-postgwas/)
[![Docs status](https://img.shields.io/github/actions/workflow/status/hlnicholls/nf-postgwas/quarto-gh-pages.yml?branch=main)](https://github.com/hlnicholls/nf-postgwas/actions/workflows/quarto-gh-pages.yml)
[![Release](https://img.shields.io/github/v/release/hlnicholls/nf-postgwas)](https://github.com/hlnicholls/nf-postgwas/releases)
[![CI](https://img.shields.io/github/actions/workflow/status/hlnicholls/nf-postgwas/smoke-test.yml?branch=main)](https://github.com/hlnicholls/nf-postgwas/actions/workflows/smoke-test.yml)



This pipeline provides a docker image with all dependencies pre-installed to run a post-genome-wide association study analysis pipeline via Nextflow.

The `nf-postgwas` pipeline integrates variant annotation, gene prioritisation, fine mapping, enrichment, and more. It is designed for flexible, modular analysis of GWAS summary statistics and downstream biological interpretation.

![nf-postgwas workflow diagram](/docs/images/nf-postgwas.png)

## Installation and Setup

First clone this repository:
```
git clone https://github.com/hlnicholls/nf-postgwas.git
```

Docker needs to be installed on your system. Instructions for installing Docker can be found [here](https://docs.docker.com/get-docker/).

**Several databases need to be mounted to the docker container to run the pipeline. These can be be downloaded from: [huggingface.co/datasets/hlnicholls/postgwas-db-grch38-2025-11](https://huggingface.co/datasets/hlnicholls/postgwas-db-grch38-2025-11)**

**Only an example 1000 genomes reference panel is given. You will also need to mount your own LD reference panel for accurate LD calculations.**

To download the huggingface dataset, you can use either Git LFS or the Hugging Face CLI/Python API.

### Download postgwas-db with Git LFS
```bash
git lfs install --skip-repo
git clone https://huggingface.co/datasets/hlnicholls/postgwas-db-grch38-2025-11
```

### Or in Python:
```python
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="hlnicholls/postgwas-db-grch38-2025-11",
    repo_type="dataset",
    local_dir="./postgwas-db-grch38-2025-11"
)
```

### Setup and run docker container

Pull the docker image and run it with the following commands (volume mounts will need to modified for your input data):

```bash
# 1) Pull the image
docker pull --platform linux/amd64 hlnicholls/postgwas-pipeline:latest

# 2) Run it with huggingface dataset downloaded and mounted
docker run -it --rm --name postgwas-pipeline \
  -v "/path/on/host/postgwas-db-grch38-2025-11:/nf-postgwas/postgwas-db-grch38-2025-11:ro" \
  -v "/path/on/host/Regenie_GWAS_Input:/nf-postgwas/GWAS_Input:ro" \
  -v "/path/on/host/postgwas_results:/nf-postgwas/results" \
  hlnicholls/postgwas-pipeline:latest
```

Conda environments in Docker image:
- ```postgwas_py39``` - primary python (v3.9) environment for most scripts. Includes PLINK2.
- ```postgwas_py27``` - python (v2.7) environment for LDSC and DEPICT.

Installed Software
- Miniconda3
- Nextflow
- nf-test
- PLINK1.9 and PLINK2
- MAGMA
- R and required R packages

## Running the Pipeline

1. Check the `/conf/local.config` memory and CPU settings (or whichever config you plan to use)
2. Set your input params in `params.yaml` in `/nf-postgwas/params.yaml` 
  - You will need to set the paths to your mounted input GWAS files and reference databases
  - Downloaded databases need to be mounted and paths set correctly
  - For several pipeline steps, you will need to also mount your LD reference panel separately (only a synthetic panel is provided in the docker image for testing)

3. Run the pipeline with Nextflow:

```bash
conda activate postgwas_py39 # needs to be activated at start
nextflow run . -profile local -params-file params.yaml
```

## Params
The main parameters to set in `params.yaml` are:
- `gwas_input_dir`: Directory containing per-chromosome GWAS files (in REGENIE format)
  - The file name must include `"_regenie_allchr.txt"` as the suffix
  - `merged_regenie_files`: Merged GWAS summary statistics file (if you have this already prepared, it must be GRCh38 and have the columns in REGENIE format)
- `output_path`: Directory for output files
- `databases`: Directory for all database subfolders
- `ld_reference_panel`: Directory for LD reference panel
- `traits`: List of trait names
- `sample_sizes`: List of sample sizes per trait

Each step in the pipeline (e.g. variant annotation, fine mapping, colocalisation) can be enabled or disabled via flags in `params.yaml`.


## Features & Workflow Steps
- **Preprocessing (PREP_FLOW):** Merges GWAS results, performs QC, liftover, and RSID mapping.
  - Merge regenie output files per chromosome into one file
    - This step is skipped if a merged GWAS summary statistics file is provided in `params.yaml` (`merged_regenie_files`).
  - Quality control filtering (MAF > 0.01 and INFO > 0.3)
  - Liftover from GRCh38 to GRCh37 - providing a "GENPOS_hg19" column in the output (downstream analyses are based on GRCh38 coordinates where possible)
  - RSID mapping using 1000 genomes phase 3 reference panel
- **Loci Identification (LOCI_FLOW):** Compiles significant loci from GWAS results.
  - Identifies genome-wide significant variants (p < 5e-8)
  - Performs distance-based pruning of significant variants to define loci (+/-500kb window)
  - Creates locus-specific GWAS summary statistics files
- **LD Annotation (LD_FLOW):** Annotates loci with linkage disequilibrium information.
  - Calculates LD metrics (r2) using user input reference panel
  - Uses PLINK to compute LD metrics:
  ```
        --ld-window-kb 4000 \
        --ld-window-r2 0.1 \
        --ld-window 99999 \
  ```
  - Annotates loci and their nearby proxies (r2 > 0.1 and +/- 500kb) with nearby genes (within 10kb)
- **Loci Grouping (LOCI_GROUPS_FLOW):** Groups loci.
  - Groups loci based on proximity (within +/-500kb) and LD (r2 > 0.4 between lead variants)
- **Variant Annotation (VARIANT_ANNOT_FLOW):** Annotates variants using CADD, pCHiC, RegulomeDB.
  - Annotates lead variants and their nearby proxies (r2 > 0.1 and +/- 500kb) with functional scores and regulatory information.
  - VEP annotation is currently not implemented in the pipeline but can be added.
  - A VCF file is output for use in the web applications of VEP or CADD.
- **MAGMA Gene Analysis (MAGMA_FLOW):** Performs gene-based association analysis.
- **PoPS Prioritisation (POPS_FLOW):** Prioritises genes using PoPS, requires MAGMA results.
- **Plotting (PLOT_FLOW):** Generates summary plots for GWAS and loci.
  - Manhattan, QQ, and QC plots for GWAS results.
  - LocusZoom plots for each identified locus.
- **LDSC (LDSC_FLOW):** Runs LD Score Regression for heritability and genetic correlation.
  - Requires >1 trait to run genetic correlation.
  - Uses precomputed LD scores (UKBB included in `/databases/ldsc` - this needs to be changed as required by the user).
- **Colocalisation (COLOC_FLOW):** Performs colocalisation analyses, including GTEx, and a custom input trait.
  - User provides a custom trait GWAS summary statistics file for colocalisation with input GWAS loci.
  - GTEx v8 eQTL summary statistics need to be included in `/databases/coloc_gwas/gtex_v8_eqtl/` (not included for all tissues in the docker image due to size).
- **Fine Mapping (FINE_MAPPING_FLOW):** Fine-maps loci to pinpoint causal variants (method: wakefield - https://annahutch.github.io/corrcoverage/, other methods in development).
- **Enrichment (ENRICHMENT_FLOW):** Runs gene set enrichment (enrichR).
- **Prioritisation (PRIORITISATION_FLOW):** Integrates results for gene prioritisation.
  - A simple weighting criteria (between 1-4) is applied to prioritise genes within loci based on results from PoPS, Coloc, HiC, GTEx (enrichR), relevant Mendelian disease genes in OMIM (enrichR), or ClinGen, and relevant IMPC mouse phenotypes. Relevant disease terms need to be provided in `params.yaml`.
  
  Each annotation is weighted equally, with (1) OMIM and ClinGen considered jointly (i.e., the presence of either counts as one scored annotation) and (2) mouse models, (3) GTEx and Hi-C also considered jointly in the same way. Finally, (4) if a gene is in a colocalising locus from `COLOC_FLOW` this will also be annotated. Higher scores indicate stronger multi-source support for gene prioritisation. For example, a score of 4 is the highest possible score (a gene has all four lines of evidence), and a score of one is the lowest (a gene only has one line of supporting evidence).
  
- **Druggability (DRUGGABILITY_FLOW):** Assesses druggability of prioritised genes.
  - Uses OpenTargets to annotate gene-drug targets. Identifies relevant drugs for prioritised genes via input disease terms of interest in `params.yaml`.
  - Uses DGIdb data to annotate identify potentially "druggable" gene targets.

## Documentation (GitHub Pages)

Full documentation is published to GitHub Pages and is available at:

https://hlnicholls.github.io/nf-postgwas/
