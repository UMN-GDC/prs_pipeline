# sandbox_singularity_runner.sh

The **primary wrapper** for running the Single Ancestry PRS Pipeline inside a **Singularity container** (`prsv2_latest.sif`). Designed for both SLURM HPC and local execution, with automatic bind-mount detection that only mounts directories present on the host system.

## Overview

The script:
1. Parses command-line flags and a shell-readable config file
2. Auto-discovers filesystem paths from the config and binds their top-level directories into the container
3. Launches `run_single_ancestry_PRS_pipeline.sh` inside the SIF via `singularity exec`

All dependencies (R 4.3.0, bigsnpr, PLINK, PRSice-2, etc.) live inside the container — nothing needs to be installed on the host except Singularity/Apptainer.

**PRSice-2 test evaluation** differs from C+T/LDpred2/lassosum2: the pipeline applies training-derived parameters (SNP set from PLINK clumping + best p-value threshold from PRSice2 `.summary`) to test data via PLINK `--score` + `--extract`, rather than re-running PRSice2 from scratch. This mirrors the C+T evaluation strategy and prevents model selection leak from test data.

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

| Variable | Required | Description |
|----------|----------|-------------|
| `summary_stats_file` | yes | GWAS summary statistics (see [Input File Formats](#input-file-formats) below) |
| `summary_stats_files` | for multi | Comma-separated list for multi-phenotype mode (overrides `summary_stats_file`). Single entry + `multi_pheno_file` also triggers multi-mode |
| `multi_pheno_file` | for multi | Multi-phenotype file with header (see [Input File Formats](#input-file-formats)) |
| `phenotype_info_file` | no | External phenotype file (see [Input File Formats](#input-file-formats)); header **auto-detected and prepended if missing** (checks if first field is `FID`/`fid`) |
| `bim_file_path` | yes | `.bim` file for allele alignment |
| `study_sample` | yes | PLINK prefix (no extension) for the study cohort |
| `output_path` | yes | Directory for all results |
| `path_repo` | yes | Path to the cloned `prs_pipeline` repository |
| `gwas_pca_eigenvec_file` | no | PCA eigenvector file for covariates |
| `afreq_file` | **required** for LDpred2/lassosum2 | PLINK2 `.afreq` file; avoids unreliable MAF computation from genotype matrix |
| `n_total_gwas` | no (default: 31968) | GWAS sample size (used for `n_eff`) |
| `skip_ss_generation` | no (default: 0) | Set `1` to skip `prepare_sumstats.R` alignment |
| `ncores` | no (default: 16) | CPU cores for parallel LD computation. Must match Slurm `--cpus-per-task` |
| `ld_cache_dir` | no | Per-chromosome LD matrix cache directory (avoids recomputing on re-runs) |
| `ld_matrix_dir` | no | Pre-computed LD matrix directory (from `src/generate_ld_matrix.R`); takes priority over `ld_cache_dir` |
| `binary_flag` | no (default: F) | Set `T` for binary (case/control) phenotypes, `F` for quantitative |
| `test_sample` | no | PLINK prefix for held-out test data (triggers `score_test.sh` after training) |
| `test_pca_eigenvec_file` | no | PCA eigenvector file for test evaluation; if omitted, computed via `plink --pca 6` |
| `RUN_CT` | no | Set `true` to run C+T |
| `RUN_LDPRED2` | no | Set `true` to run LDpred2 |
| `RUN_LASSOSUM2` | no | Set `true` to run lassosum2 |
| `RUN_PRSice2` | no | Set `true` to run PRSice-2 (requires C+T) |

### Config-based method control

When `RUN_*` variables are defined in the config, no CLI method flags (`-c -l -s -P`) are needed — the config values take precedence (sourced before arg parsing). Example:

```bash
# Inside config.txt — controls which methods run
RUN_CT=true
RUN_LDPRED2=true
RUN_LASSOSUM2=false
RUN_PRSice2=true
```

## Input File Formats

All input files are tab-separated (`.tsv`). Headers are **required** (auto-detected and prepended for `phenotype_info_file`; see below).

### Summary Statistics (`summary_stats_file`)

GWAS summary statistics with at minimum these columns:

```text
SNP	CHR	BP	A1	A2	beta	beta_se	P	n_eff
rs101	1	10001	G	A	0.012	0.005	0.0081	50000
rs102	1	10005	C	T	-0.008	0.004	0.0234	50000
```

Column mapping is flexible: `rsid`/`rs_id`/`rsids` for SNP, `pval` for P. The first step (`prepare_sumstats.R`) expects `beta` (lowercase) and `sebeta`/`beta_se` for SE. When `skip_ss_generation=1`, downstream R scripts also accept `BETA`/`eff` and `SE`.

### External Phenotype File (`phenotype_info_file`)

Single-phenotype file with 3+ columns. The pipeline uses **column 3** (first phenotype column) for analysis.

```text
FID	IID	phenotype
Sample1	ID001	1.24
Sample2	ID002	-0.87
Sample3	ID003	0.53
```

**Feature: Auto-Header Detection**
- If the file lacks a header (first field ≠ `FID`/`fid`), the pipeline auto-detects the column count and prepends `FID	IID	phenotype1	[phenotype2]	...`
- Works for **any number of columns** (1 phenotype or many)
- A log message confirms: `No header detected — prepending N-column header (FID, IID, phenotype1...phenotypeN)`
- **No requirement change:** the single-phenotype pipeline only uses column 3 regardless of how many columns exist. Use `multi_pheno_file` (see below) to run multiple phenotypes in a single pass.

### Multi-Phenotype File (`multi_pheno_file`)

Tab-separated file with FID, IID, and one or more phenotype columns. **Header is required** for column name lookup.

```text
FID	IID	phenoA	phenoB	phenoC
Sample1	ID001	1.24	0.53	-0.10
Sample2	ID002	-0.87	1.12	0.75
```

The pipeline extracts each named column into individual `study_sample_pheno.txt` files and runs all requested methods per phenotype.

**Trigger conditions:** multi-phenotype mode activates when `summary_stats_files` has ≥2 comma-separated entries OR a single entry combined with `multi_pheno_file`.

**1:1 matching required:** the number of summary stats files must equal the number of phenotype value columns (columns 3+). For example, 3 phenotype columns require 3 summary stats files:

```bash
summary_stats_files="/path/to/gwas_A.txt,/path/to/gwas_B.txt,/path/to/gwas_C.txt"
multi_pheno_file="/path/to/phenotypes.txt"
```

The pipeline validates this at runtime and exits with an error if the counts differ. Output goes to `${output_path}/prs_pipeline/<PHENO_NAME>/<METHOD>/`.

### Study Sample PLINK Files (`study_sample`)

Standard PLINK binary format (`.bed`/`.bim`/`.fam`). The `.fam` file provides the baseline phenotype (column 6) — used as a fallback when no external phenotype file is given.

### PCA Eigenvectors (`gwas_pca_eigenvec_file`)

PLINK2-format eigenvector file — 2 header columns (FID, IID) followed by PCs:

```text
FID	IID	PC1	PC2	PC3	PC4	PC5	PC6
Sample1	ID001	-0.021	0.015	0.003	-0.008	0.011	-0.004
```

### Allele Frequencies (`afreq_file`)

PLINK2 `.afreq` format:

```text
CHR	SNP	REF	ALT	ALT_FREQS	OBS_CT
1	rs101	G	A	0.350	50000
1	rs102	C	T	0.620	50000
```

**Required** for LDpred2 and lassosum2 to avoid unreliable MAF computation from the genotype matrix.

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
    CT_PRSice2_summary_stat_file.txt    # Aligned summary stats (after prepare_sumstats.R)
    study_sample_pheno.txt              # Phenotype file
  prs_pipeline/
    CT/              -- Clumping + Thresholding results
    LDpred2/         -- Bayesian PRS (infinitesimal + grid models)
    lassosum2/       -- Penalized regression grid results
    PRSice2/         -- PRSice-2 outputs (best, summary, plots)
    logs/            -- PRSice2.log (redirected here; last 50 lines on failure)
    test_evaluation/ -- Present only if test_sample was configured
                       <method>_results.txt, <method>_scores.txt
```

**Multi-phenotype** adds a phenotype-name subdirectory: `${output_path}/prs_pipeline/<PHENO_NAME>/<METHOD>/`.

**CPU note:** the sandbox runner Slurm header requests **4 CPUs** while the pipeline defaults to `ncores=16`. Set `ncores` in your config to match `--cpus-per-task` (e.g., `ncores=4`) to avoid CPU oversubscription.

## Sandbox vs. Production Runner

`singularity_runner.sh` is **deprecated** — use `sandbox_singularity_runner.sh` for all new work.

| Aspect | `sandbox_singularity_runner.sh` **(primary)** | `singularity_runner.sh` (deprecated) |
|--------|-----------------------------------------------|--------------------------------------|
| Container SIF | `prsv2_latest.sif` (pull: `apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest`) | `singleprs_latest.sif` |
| Bind strategy | Top-level roots, existence-checked | Mount-point extraction, always binds `/projects` |
| SLURM header | Present (20h, 4c, 64g, `medium` partition) | Absent (manual/local use only) |
| Use case | Safer for shared/transient filesystems | Legacy — hardcoded for specific HPC environment |
