# prs_pipeline
A repository to aid in the use of multiple PRS models. 

Assumes the inputs are plink formatted and split by ancestry. Also assumes they are providing pc calculations as an input

## Helper resource
Sequence to run scripts if using simulated data and needing to generate summary stats files
1. split_top_n_subjs.sh
2. run_split_plink_data.sh
3. generate_summary_stat_files.sh
4. restructure_output_dir.sh

run_prepare_prs.sh script runs the four modules to prepare for running a prs method in sequence for you.

```{bash}
Sample call with a smaller data set
sbatch prs_pipeline/run_prepare_prs.sh -1 /scratch.global/baron063/simulations/temp/AFR_simulation -2 /scratch.global/baron063/simulations/temp/EUR_simulation -P /projects/standard/gdc/public/prs_methods/data/adjusted_1kgPCs.tsv -n 300
```

## PRS-CSx Joint Ancestry Pipeline
This script runs a joint-ancestry PRS analysis using PRS-CSx, starting from GWAS summary statistics and PLINK genotype files, and ending with ancestry-specific PRS scores and R squared evaluation. It is designed to run on an HPC cluster using Slurm.

The original toolset can be found here [PRS-CSx](https://github.com/getian107/PRScsx). 

*The development and evaluation of PRS-CSx are described in:
Y Ruan, YF Lin, YCA Feng, CY Chen, M Lam, Z Guo, Stanley Global Asia Initiatives, L He, A Sawa, AR Martin, S Qin, H Huang, T Ge. Improving polygenic prediction in ancestrally diverse populations. Nature Genetics, 54:573-580, 2022.*

### Overview
At a high level the pipeline:
1. Validates inputs and optional configuration files
2. Reformats GWAS summary statistics into PRScsx-compatible format
3. Runs `PRScsx.py` jointly across two ancestries
4. Combined chromosome-level PRScsx outputs into genome-wide SNP weight files
5. Computes PRS scores using `plink2`
6. Evaluates predictive performance (R squared and adjusted R squared) using an R helper script

### Requirements
* Bash (with Slurm support)
* Python (compatible with PRS-CSx)
* PRS-CSx GitHub repository cloned [PRS-CSx](https://github.com/getian107/PRScsx)
* PRS-CSx LD reference panels downloaded [PRS-CSx LD panels](https://github.com/getian107/PRScsx/blob/master/README.md)
* Plink2 [plink2.0 home](https://www.cog-genomics.org/plink/2.0/)
* R (for PRS evaluation)
* Conda environment loaded via:
```bash
source /projects/standard/gdc/public/envs/load_miniconda3.sh
```

### Usage
```bash
sbatch run_PRScsx.sh \
  --path_code PATH \
  --path_data_root PATH \
  --path_ref_dir PATH \
  --path_plink2 PATH \
  --anc1 AFR \
  --anc2 EUR \
  --target_sumstats_file PATH \
  --training_sumstats_file PATH \
  --output_dir PATH \
  --reference_SNPS_bim PATH \
  --study_sample_plink PATH \
  --study_sample_plink_anc2 PATH
```

Alternatively, you may supply a shell-readable configuration file:
```bash
sbatch run_PRScsx.sh --c config.sh
```
When `--c` is provided, *all defaults and command-line arguments are overridden* by values defined in the config file.

### Key Arguments
| Argument                    | Description                                       |
| --------------------------- | ------------------------------------------------- |
| `--path_code`               | Path to PRScsx source code (contains `PRScsx.py`) |
| `--path_data_root`          | Root directory for GWAS and genotype data         |
| `--path_ref_dir`            | PRScsx reference LD directory (e.g. 1KG)          |
| `--path_plink2`             | Path to `plink2` binary                           |
| `--anc1`                    | Target ancestry (e.g. AFR)                        |
| `--anc2`                    | Training/reference ancestry (e.g. EUR)            |
| `--target_sumstats_file`    | GWAS summary stats for target ancestry            |
| `--training_sumstats_file`  | GWAS summary stats for training ancestry          |
| `--reference_SNPS_bim`      | BIM prefix for reference SNP set                  |
| `--study_sample_plink`      | PLINK prefix for ancestry 1 samples               |
| `--study_sample_plink_anc2` | PLINK prefix for ancestry 2 samples               |
| `--output_dir`              | Base directory for all outputs                    |

### Input Summary Statistics Format
The input GWAS summary statistics are expected to include:
* SNP ID
* Effect allele (A1)
* Non-effect allele (A2)
* Effect size (BETA)
* Standard error (SE)
* Sample size (used to extract `n_gwas`)
The script reformats these into PRScsx-compatible files:

### Outputs
All outputs are written to:
```bash
${output_dir}/prs_pipeline/PRScsx/
```
Key output files include: 
* Combined SNP weights
    * `PRScsx_<ANC>_combined_weights.txt`
* PRS scores
    * `PRScsx_joint_<ANC>_score.sscore`
* Model performance
    * `<ANC>_PRS_sscore_Rsqr.txt`
    * `<ANC>_adj_PRS_sscore_Rsqr.txt`

Each ancestry is processed separately after the joint PRScsx run, but shares the same fitted model.

### Notes
* The script assumes autosomes only (chromosomes 1-22).
* Random seed is fixed (`--seed=42`) for reproducibility.
* Designed for simulated or real multi-ancestry GWAS with matched SNP sets.

### Sample config file
```bash
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx" # Path to cloned PRS-CSx GitHub repository
path_data_root="/projects/standard/gdc/public/prs_methods/data/simulated_1000G" # Path to genomic data 
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref" # Path to the downloaded reference LD panels 
path_plink2="/projects/standard/gdc/public/plink2" # Path to plink2 executable
anc1="AFR" # Ancestry 1 # Must be one of the five super populations
anc2="EUR" # Ancestry 2
target_sumstats_file="${path_data_root}/gwas/target_sumstats.txt" # Summary stats file for Anc1
training_sumstats_file="${path_data_root}/gwas/training_sumstats.txt" # Summary stats file for Anc2
output_dir="${path_data_root}" # Desired place to output to
reference_SNPS_bim="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample" # Full path to the reference SNPs without the extension # Assumed that there is a .bim file
study_sample_plink="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample" # Full path to the study population plink files for Ancestry 1
study_sample_plink_anc2="${path_data_root}/anc2_plink_files/${anc2}_simulation_study_sample" # Full path to the study population plink files for Ancestry 2
prs_pipeline="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline" # Full path to the cloned prs_pipeline GitHub repository
```

## VIPRS PRS Pipeline

This script runs a complete **VIPRS-based polygenic risk score (PRS) pipeline** on an HPC system using Slurm. It generates GWAS summary statistics from PLINK genotype data, fits a VIPRS model, computes PRS scores in a study sample, and evaluates predictive performance.

The original toolset can be found [Variational Inference of Polygenic Risk Scores](https://shz9.github.io/viprs/).

*Zabad, S., Gravel, S., & Li, Y. (2023). Fast and accurate Bayesian polygenic risk modeling with variational inference. The American Journal of Human Genetics, 110(5), 741-761. https://doi.org/10.1016/j.ajhg.2023.03.009*


### What this script does

At a high level, the script:

1. Activates a Conda environment containing VIPRS
2. Runs `plink2 --glm` by chromosome to generate GWAS summary statistics
3. Reformats PLINK output into VIPRS-compatible summary statistics
4. Fits a VIPRS model using LD reference panels
5. Computes PRS scores in an independent study sample
6. Evaluates PRS performance using phenotype and covariate files

### Requirements

- Slurm-enabled HPC environment
- Bash
- PLINK 2
- VIPRS (`viprs_fit`, `viprs_score`, `viprs_evaluate`)
- Conda environment containing VIPRS
- LD reference panels for VIPRS

The script loads the environment using:
```bash
source /projects/standard/gdc/public/envs/load_miniconda3.sh
conda activate viprs_env
```

### Usage
This script requires a configuration file. Run using `sbatch`:
```bash
sbatch run_viprs.sh --c config.sh
```
The configuration file is a shell-readable file that defines input paths, covariates, and output locations. All defaults in the script are overridden by values in the config file.

#### Sample configuration file

```bash
path_data=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
out_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_plink2=/projects/standard/gdc/public/plink2
bfile_gwas_input=${path_data}/anc1_plink_files/archived/AFR_simulation_gwas
bfile_study_sample=${path_data}/anc1_plink_files/AFR_simulation_study_sample
covariate_file_gwas=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/prs_pipeline/viprs/gwas/temp/viprs_summary_stats_covar_sex_no_header.txt
covariate_file_study_sample=${path_data}/prs_pipeline/viprs/study_sample_covar.txt
```

### Inputs
* PLINK genotype files for GWAS (--bfile format)
* PLINK genotype files for the study/sample cohort
* Covariate file for GWAS (used in plink2 --glm)
* Covariate file for PRS evaluation
* VIPRS LD reference panels (split by chromosome) downloadable following the original documentation [Download LD Reference](https://shz9.github.io/viprs/download_ld/)

### Outputs
All outputs are written to the VIPRS output directory:
* GWAS summary statistics (per chromosome)
* VIPRS posterior effect size estimates:
    * `VIPRS_EM.fit.gz`
* Polygenic scores:
    * `VIPRS_PGS.prs`
* Evaluation results:
    * `viprs_evaluate_results.eval`
Intermediate PLINK logs and temporary files are organized into subdirectories automatically.

### Notes
* Analysis is restricted to autosomes (chromosomes 1-22)
* GWAS summary statistics are generated using an additive genetic model
* Phenotype files are generated automatically from provided study population plink files. It is possible to adjust this if desired.
* Designed for simulated or real genotype data split by ancestry

## Single Ancestry PRS Pipeline
A modular SLURM-based pipeline for calculating Polygenic Risk Scores (PRS) using multiple state-of-the-art methods. This tool automates the heavy lifting of data alignment and allows you to run several PRS algortihms in parallel.

### Features
* **Multi-Method Support**: Run **C+T**, **LDpred2**, **lassosum2**, and **PRSice-2** with a single command.
* **Automated Alignment**: Handles the synchronization between GWAS summary statistics and target genotype `.bim` files.
* **HPC Optimized**: Build-in SLURM directives and parallel background processing for efficient resource usage.
* **Configurable**: Supports external configuration files to keep your scripts clean and reporducible.

### Prerequisites
**HPC Evnironment**
* Scheduler: SLURM
* R Version: 4.4.0+ (OpenBLAS recommended)
* R Packages: `bigsnpr `, `optparse`, `data.table`, `magrittr  `

**Required Tools**
* PRSice-2: The executable should be located within your `${path_repo}/src/` directory.

**Summary Statistics File Format**
The file needs to have a column header and contain the following columns. If there are multiple listed in `{}` that means that any 1 of the options is recognized but the file cannot contain more than 1 of these overlapping columns.
* {rsid, rs_id, rsids}
* {A1, alt, a1}
* {p, pval}
* {sebeta, beta_se}
* beta
* ref
* chrom

### Usage

#### Basic Execution
Submit the job to the SLURM queue by selecting your desired methods via flags including the path to your config file.

```bash
sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh -c -l -s -P -C config.txt
```

#### Flag Reference
| Flag | Description |
| -------- | ------- |
|  `-C <file>` | **Config**: load an external file to override default paths. |
| `-c` | **C+T**: Run Clumping + Thresholding. |
| `-l` | **LDpred2**: Run the LDpred2 algorithm (via `bigsnpr`) |
| `-s` | **lassosum2**: Run the lassosum2 algorithm (via `bigsnpr`) |
| `-P` | **PRSice2**: Run PRSice-2. Note: Requires `-c` to be active. |
| `-B` | **Binary**: Use this flag for binary phenotypes (Case/Control). |
| `-S` | **Skip**: Skip the initial summary stats alignment step (Use with caution). |

### Configuration
Instead of modifying the main script, create a `my_project.conf` file to define your data paths:

```bash
# my_project.conf
summary_stats_file="/path/to/your/gwas_stats.tsv"
bim_file_path="/path/to/your/genotypes.bim"
study_sample="/path/to/your/genotype_prefix"
output_path="/path/to/results"
n_total_gwas=31968
```
Then run:
```
sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh -c -l -s -P -C config.txt
```

### Output Structure
The pipeline automatically organizes results into a clean directory tree:

```Plaintext
output_path/
|-- gwas/
|   --- CT_PRSice2_summary_stat_file.txt  # Aligned stats
|   --- study_sample_pheno.txt            # Extracted phenotypes
|-- prs_pipeline/
|   --- CT/          # Clumping + Thresholding results
|   --- LDpred2/     # Bayesian PRS results
|   --- lassosum2/   # Penalized regression results
|   --- PRSice2/     # PRSice-2 tables and plots
```

### Imporant Notes
* **Resource Allocation**: The script defaults to **16 CPUs** and **64GB RAM**. Adjust the `#SBATCH` headers if your LD reference panel or genotype file is exceptionally large.
* **C+T Dependency**: The PRSice-2 implementation in this script relies on the data preparation steps performed during the C+T run. Always include `-c` when using `-P`.
* **Environment**: Ensure `R_LIBS_USER` in the script points to the library where `bigsnpr` is installed.

### Troubleshooting
Memory Errors: If LDpred2 fails, ensure the --mem=64g SLURM header is sufficient for your LD reference panel.
Missing Variants: If the aligned summary stats file is empty, check that your .bim file RSIDs match the format in your summary statistics.

### Original References and Documentation
The following methods are implemented in this pipeline. If you use these results in a publication, please cite the corresponding papers:
1. Clumping + Thresholding (C+T)
* Original Method: The foundational approach used since the early days of GWAS (e.g., Purcell et al., 2009).
* Key Reference: [Choi et al. (2020) - A guide to performing Polygenic Risk Score analyses](https://doi.org/10.1038/s41596-020-0353-1)
* Implementation Note: This pipeline uses a standard clumping approach typically executed via PLINK.
2. LDpred2
* Original Paper: [Privé et al. (2020) - LDpred2: better, faster, stronger](https://doi.org/10.1093/bioinformatics/btaa1029)
* Documentation: [bigsnpr - LDpred2 Tutorial](https://privefl.github.io/bigsnpr/articles/LDpred2.html)
* GitHub: [privefl/bigsnpr](https://github.com/privefl/bigsnpr)
3. lassosum2
* Original Paper: [Privé et al. (2022) - lassosum2: an updated version complementing LDpred2](https://www.google.com/search?q=https://doi.org/10.1016/j.xhgg.2022.100136)
* Documentation: [bigsnpr - lassosum2 Reference](https://privefl.github.io/bigsnpr/reference/snp_lassosum2.html)
* GitHub: [privefl/bigsnpr](https://github.com/privefl/bigsnpr)
4. PRSice-2
* Original Paper: [Choi and O'Reilly (2019) - PRSice-2: Polygenic Risk Score software for biobank-scale data](https://doi.org/10.1093/gigascience/giz082)
* Documentation: [Official PRSice-2 Documentation](https://www.prsice.info/)
* GitHub: [ChoiS_github/PRSice](https://github.com/choishingwan/PRSice)

#### Quick Comparison of Methods

| Method | Approach | LD Handling |
| -------- | ------- | ------- |
| C+T | Heuristic | Physical/Correlation Clumping |
| LDpred2 | Bayesian | Gibbs Sampler (Markov Chain Monte Carlo) |
| lassosum2 | Penalized Regression | Elastic Net / Coordinate Descent |
| PRSice-2 | C+T Optimization | Automated High-Resolution Clumping |

## TL-PRS
* TL-PRS - plink files are separate 

## CTSLEB
* CTSLEB - plink files are together for anc1 and anc2, They are separated for summary statitics and reference panels. 

## Prosper
* prosper - plink files are together for anc1 and anc2

# Single ancestry PRS modules
