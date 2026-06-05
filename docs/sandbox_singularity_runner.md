# sandbox_singularity_runner.sh

The **primary wrapper** for running the Single Ancestry PRS Pipeline inside a **Singularity container** (`prsv2_latest.sif`). Designed for both SLURM HPC and local execution, with automatic bind-mount detection that only mounts directories present on the host system.

## Overview

The script:
1. Parses command-line flags and a shell-readable config file
2. Auto-discovers filesystem paths from the config and binds their top-level directories into the container
3. Launches `run_single_ancestry_PRS_pipeline.sh` inside the SIF via `singularity exec`

All dependencies (R 4.3.0, bigsnpr, PLINK, PRSice-2, etc.) live inside the container — nothing needs to be installed on the host except Singularity/Apptainer.

## Prerequisites

- **Singularity** (or Apptainer) installed on the host
- **`prsv2_latest.sif`** — the pre-built container image, download instructions below
- A **config file** (see [Configuration](#configuration)) pointing to existing paths on the host

### Download Singularity Image Instructions

```bash
apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest
```

## Usage

```bash
# SLURM submission
sbatch sandbox_singularity_runner.sh --C config.txt -c -l -s -P

# Local execution
bash sandbox_singularity_runner.sh --C config.txt -c -l
```

### Flags

| Flag | Description |
|------|-------------|
| `-C <file>` / `--C <file>` | Path to shell-readable config file (required) |
| `--bind <paths>` | Comma-separated bind paths (overrides auto-detection) |
| `-c` | Run Clumping + Thresholding (C+T) |
| `-l` | Run LDpred2 |
| `-s` | Run lassosum2 |
| `-P` | Run PRSice-2 *(requires `-c`)* |
| `-S` | Skip summary statistics alignment (use with caution) |
| `-B` | Binary phenotype (case/control); omit for quantitative traits |

## Configuration

The config file is a shell-readable file sourcing variables directly into the script. A full template is at `templates/single_anc_config.txt`.

```bash
summary_stats_file="/path/to/gwas_stats.tsv"
bim_file_path="/path/to/genotypes.bim"
study_sample="/path/to/genotype_prefix"
output_path="/path/to/results"
path_repo="/path/to/prs_pipeline"
gwas_pca_eigenvec_file="/path/to/pca.eigenvec"
afreq_file="/path/to/sample.afreq"
```

| Variable | Description |
|----------|-------------|
| `summary_stats_file` | GWAS summary statistics (see README for required columns) |
| `summary_stats_files` | Comma-separated list for multi-phenotype mode (overrides `summary_stats_file`) |
| `multi_pheno_file` | Tab-sep file: FID, IID, pheno1, pheno2, ... (with header) for multi-phenotype mode |
| `phenotype_info_file` | External phenotype file (FID IID Value, no header) |
| `bim_file_path` | `.bim` file for allele alignment |
| `study_sample` | PLINK prefix (no extension) for the study cohort |
| `output_path` | Directory for all results |
| `path_repo` | Path to the cloned `prs_pipeline` repository |
| `gwas_pca_eigenvec_file` | PCA eigenvector file for covariates |
| `afreq_file` | PLINK2 `.afreq` file; **required** for LDpred2/lassosum2 to avoid unreliable MAF computation from genotype matrix |
| `ld_cache_dir` | Per-chromosome LD matrix cache directory (avoids recomputing on re-runs) |
| `ld_matrix_dir` | Pre-computed LD matrix directory (from `src/generate_ld_matrix.R`) |

## Auto-Bind Mount Detection

The sandbox runner reads every config path variable and:

1. Strips quotes and extracts the path value
2. Extracts the **top-level root directory** (e.g., `/scratch.global` from `/scratch.global/baron063/...`)
3. Only binds that root if it **actually exists** on the host (`[[ -d "$root_dir" ]]`)
4. Deduplicates so each root is bound only once

This is more cautious than the production runner — it never binds a directory that isn't present on the current machine. Non-existent paths produce a warning to stderr and are silently skipped.

## Outputs

Results are written under `output_path/`:

```
output_path/
  gwas/
    CT_PRSice2_summary_stat_file.txt
    study_sample_pheno.txt
  prs_pipeline/
    CT/              -- Clumping + Thresholding
    LDpred2/         -- Bayesian PRS (infinitesimal + grid models)
    lassosum2/       -- Penalized regression grid results
    PRSice2/         -- PRSice-2 outputs (best, summary, plots)
```

## Sandbox vs. Production Runner

`singularity_runner.sh` is **deprecated** — use `sandbox_singularity_runner.sh` for all new work.

| Aspect | `sandbox_singularity_runner.sh` **(primary)** | `singularity_runner.sh` (deprecated) |
|--------|-----------------------------------------------|--------------------------------------|
| Container SIF | `prsv2_latest.sif` (pull: `apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest`) | `singleprs_latest.sif` |
| Bind strategy | Top-level roots, existence-checked | Mount-point extraction, always binds `/projects` |
| SLURM header | Present (20h, 16c, 200g, `medium` partition) | Absent (manual/local use only) |
| Use case | Safer for shared/transient filesystems | Legacy — hardcoded for specific HPC environment |
