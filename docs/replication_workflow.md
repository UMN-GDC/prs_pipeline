# End-to-End PRS Replication Workflow

This document describes how to replicate a single-ancestry PRS analysis starting from raw ABCD/UK Biobank data, using the `prs_pipeline` toolchain on MSI HPC with Singularity/Apptainer.

## Overview

```
raw genotype data  ──>  genomic_preps.sh  ──>  cleaned PLINK files + .afreq
raw sumstats       ──>  transform columns  ──>  cleaned_sumstats.tsv
                                                    │
                        config file (paths to all of the above)
                                    │
                                    └──>  sbatch sandbox_singularity_runner.sh --C config.txt
                                              │
                                              ├── prepare_sumstats.R (inside container)
                                              ├── C+T
                                              ├── LDpred2
                                              ├── lassosum2
                                              └── PRSice-2
```

## 1. Summary Statistics

### Source

The raw GWAS summary statistics come from the **UK Biobank BIG40** resource:

> https://open.oxcin.ox.ac.uk/ukbiobank/big40/BIG40-IDPs_v4/IDPs.html

### Required Transformations

The raw file must be cleaned before `prepare_sumstats.R` can read it. Apply these steps with `awk`:

```bash
# Strip leading zeros from chromosome column ($1), convert -log10Pval ($NF) to pval,
# rename column 7 header to beta_se
awk '
BEGIN {OFS="\t"}
NR==1 {
  $7 = "beta_se"
  $NF = "pval"
  print
  next
}
{
  gsub(/^0+/, "", $1)          # remove leading zeros: "01" -> "1"
  if ($1 == "") $1 = "0"       # if all zeros, set to "0"
  $NF = 10^(-$NF)              # convert -log10(P) to pval
  print
}' raw_sumstats.tsv > cleaned_sumstats.tsv
```

At this point the cleaned file is ready. You do **not** run `prepare_sumstats.R` yourself — it runs automatically inside the container as the first step of the pipeline when `skip_ss_generation=0` in the config. That script:

1. Reads the cleaned tab-separated sumstats
2. Detects required columns by flexible name matching (e.g. `rsid`/`rs_id`/`rsids`, `A1`/`alt`/`a1`, `p`/`pval`, `sebeta`/`beta_se`)
3. Joins with the target `.bim` file to align alleles and drop unmatched variants
4. Adds `n_eff` (from `--n_total` or estimated from allele frequency and SE)
5. Deduplicates variants and writes a 9-column file to `gwas/CT_PRSice2_summary_stat_file.txt`

**Expected output columns**: `SNP`, `CHR`, `BP`, `A1`, `A2`, `beta`, `beta_se`, `P`, `n_eff`

## 2. Phenotype File and Gender File

These two files must be prepared **before** running `genomic_preps.sh`. Both are plain-text, tab-separated, **without a header row**, with the IID duplicated in columns 1 and 2 (FID and IID are the same — there are no family groupings). **Ensure there are no quotes in the file for FID and IID**

### Gender file (`ready_sex.txt`)

Format: `FID IID <sex_code>` where sex code is `1` (male) or `2` (female) per PLINK convention.

```bash
# Sex column from a CSV: extract IID and sex column, duplicate IID as FID
awk -F, 'NR>1 {print $1, $1, $2}' cohort_sex.csv > ready_sex.txt
#                                          ↑ IID   ↑ sex
```

### Phenotype file (`ready_pheno.txt`)

Format: `FID IID <phenotype_value>`

```bash
# Phenotype from a CSV with IID in column 1 and the measure in column 3
awk -F, 'NR>1 {print $1, $1, $3}' idp_values.csv > ready_pheno.txt
#                                          ↑ IID   ↑ phenotype
```

## 3. Ancestry Subsetting

Ancestry keep-lists are plain-text files with one `FID IID` pair per line, generated from a CSV that assigns each IID to a superpopulation. Duplicate the IID as FID.

```bash
# Extract EUR individuals from a CSV: IID + superpopulation column
awk -F, 'NR>1 && $3 == "EUR" {print $1, $1}' ancestry_assignments.csv > keep_eur.txt
#                               ↑ IID (as FID and IID)

awk -F, 'NR>1 && $3 == "AFR" {print $1, $1}' ancestry_assignments.csv > keep_afr.txt
```

## 4. genomic_preps.sh

### What it does

The `genomic_preps.sh` script runs three sequential PLINK steps inside an Apptainer container:

1. **Update phenotype and sex** — attaches the phenotype and gender files to the raw PLINK dataset
2. **Split by ancestry + QC** — extracts EUR and AFR subsets, filters variants with `--geno 0.05` and samples with `--mind 0.05`
3. **Calculate allele frequencies** — runs `plink2 --freq` on each ancestry subset, producing `.afreq` files. These are **required** for LDpred2 and lassosum2 to ensure the pipeline works as intended. Providing `.afreq` avoids MAF computation directly from the genotype matrix, which can produce unreliable results on data with missingness or insufficient variation.

### Download Singularity Image Instructions

This singularity image is needed to run genomic_preps.sh.

``` bash
apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest
```

### Before running

Edit the top section of `genomic_preps.sh` to match your paths:

```bash
CONTAINER="/home/user/prs_pipeline/prsv2_latest.sif"
RAW_DATA="/path/to/merged_chroms"                # PLINK prefix (no extension)
PHENO_FILE="/path/to/ready_pheno.txt"             # FID IID Value (no header)
GENDER_FILE="/path/to/ready_sex.txt"              # FID IID SexCode (no header)
```

And ensure your ancestry keep-lists (`keep_eur.txt`, `keep_afr.txt`) exist at the paths referenced in the script.

### Outputs

| File | Description |
|------|-------------|
| `abcd_EUR_final.bed` / `.bim` / `.fam` | Cleaned EUR PLINK files |
| `abcd_AFR_final.bed` / `.bim` / `.fam` | Cleaned AFR PLINK files |
| `abcd_EUR_final.afreq` | EUR allele frequencies (for `afreq_file` config) |
| `abcd_AFR_final.afreq` | AFR allele frequencies |

### Run

```bash
# If the script has Slurm headers, use sbatch:
sbatch genomic_preps.sh

# Or run interactively if on a compute node:
bash genomic_preps.sh
```

## 5. Config File Setup

Create a shell-readable config file pointing to the outputs from `genomic_preps.sh` and the cleaned sumstats. The `RUN_*` variables in the config control which methods execute, so no CLI flags are needed.

```bash
# Paths
summary_stats_file="/path/to/cleaned_sumstats.tsv"
bim_file_path="/path/to/abcd_EUR_filtered.bim"
study_sample="/path/to/abcd_EUR_filtered"          # PLINK prefix (no extension)
output_path="/path/to/prs_results"
path_repo="/path/to/prs_pipeline"               # cloned repository
skip_ss_generation=0                             # 0 = run prepare_sumstats.R inside container
n_total_gwas=31968
ncores=16                                        # CPUs for parallel LD computation
ld_cache_dir="/path/to/ld_cache"                 # Cache per-chromosome LD matrices (avoids recomputing)

# Parameters
gwas_pca_eigenvec_file="/path/to/pca.eigenvec"  # if using PCA covariates
afreq_file="/path/to/abcd_EUR_final.afreq"       # REQUIRED for LDpred2 and lassosum2

# Method flags — these control which methods run, no -c -l -s -P flags needed on CLI
RUN_CT=true
RUN_LDPRED2=true
RUN_LASSOSUM2=true
RUN_PRSice2=true
```

### Key config variables

| Variable | Required | Description |
|----------|----------|-------------|
| `summary_stats_file` | yes | Raw (or pre-cleaned) sumstats file. If `skip_ss_generation=0`, `prepare_sumstats.R` runs inside the container to align and deduplicate |
| `bim_file_path` | yes | `.bim` matching the study sample PLINK files |
| `study_sample` | yes | PLINK prefix (no extension) for the target cohort |
| `output_path` | yes | Directory for all PRS method outputs |
| `path_repo` | yes | Path to the cloned `prs_pipeline` repository |
| `afreq_file` | **required** | `.afreq` file from `plink2 --freq`. Required for LDpred2 and lassosum2 to bypass MAF computation from the genotype matrix |
| `ncores` | no | Number of CPU cores for parallel LD computation (default: 16). Match to Slurm `--cpus-per-task` |
| `ld_cache_dir` | no | Directory for cached per-chromosome LD matrices. If the cache exists, LD loading takes seconds instead of hours. Delete and re-run to regenerate if inputs change |
| `skip_ss_generation` | no | Set to `0` (default) to let the pipeline run `prepare_sumstats.R` inside the container; set to `1` if you ran it manually |

### Method config

The `RUN_CT`, `RUN_LDPRED2`, `RUN_LASSOSUM2`, and `RUN_PRSice2` variables in the config control which methods execute. When these are defined in the config, no CLI flags (`-c -l -s -P`) are needed on the sbatch command.

## 6. Running the Pipeline

Submit via Slurm using the sandbox singularity runner:

```bash
sbatch sandbox_singularity_runner.sh --C /path/to/config.txt
```

No CLI method flags (`-c -l -s -P`) are needed — the config file's `RUN_CT`, `RUN_LDPRED2`, `RUN_LASSOSUM2`, and `RUN_PRSice2` variables control which methods execute. Only `--C` (`-C`) is required.

### Config-based method control

| Config variable | Effect |
|----------------|--------|
| `RUN_CT=true` | Run Clumping + Thresholding |
| `RUN_LDPRED2=true` | Run LDpred2 |
| `RUN_LASSOSUM2=true` | Run lassosum2 |
| `RUN_PRSice2=true` | Run PRSice-2 (requires C+T to also be enabled) |

### Output structure

```
${output_path}/
  gwas/
    CT_PRSice2_summary_stat_file.txt
    study_sample_pheno.txt
  prs_pipeline/
    CT/           Clumping + Thresholding
    LDpred2/      Bayesian PRS (inf + grid)
    lassosum2/    Penalized regression
    PRSice2/      PRSice-2 results
```

## Awk Cheat Sheet for Data Prep

```bash
# Remove header from a file
tail -n +2 file.tsv

# Print FID and IID (columns 1-2) from a FAM file
awk '{print $1, $2}' sample.fam

# Extract column 6 (phenotype) from FAM for a pheno file
awk '{print $1, $2, $6}' sample.fam > pheno.txt

# Convert -log10(P) to P
awk '{$(NF) = 10^(-$(NF)); print}' sumstats.tsv

# Strip leading zeros from chromosome
awk '{gsub(/^0+/, "", $1); print}' sumstats.tsv

# Filter to autosomes only (1-22)
awk '$1 ~ /^[0-9]+$/ && $1 >= 1 && $1 <= 22' sumstats.tsv

# Keep specific samples by ID list
awk 'NR==FNR {keep[$1$2]=1; next} keep[$1$2]' keep.txt sample.fam > subset.fam
```
