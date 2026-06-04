#!/usr/bin/env Rscript
#
# Generate test PLINK dataset + multi-phenotype inputs for PRS pipeline testing.
#
# Usage: Rscript src/generate_test_data.R <output_dir>
# Example: Rscript src/generate_test_data.R test_data
#
# Creates:
#   test_data/
#   ├── study_sample.bed       # PLINK binary (200 SNPs, 100 indiv)
#   ├── study_sample.bim       # SNP map
#   ├── study_sample.fam       # Subject info
#   ├── gwas.bim               # GWAS BIM (same as study_sample.bim)
#   ├── combined_phenotypes.txt  # FID, IID, phenoA, phenoB, phenoC
#   ├── gwas_phenoA.txt        # Summary stats for phenoA
#   ├── gwas_phenoB.txt        # Summary stats for phenoB
#   └── gwas_phenoC.txt        # Summary stats for phenoC
#   └── test_config.txt        # Ready-to-use config for run_single_ancestry_PRS_pipeline.sh

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  out_dir <- "test_data"
  message("No output dir specified; using: ", out_dir)
} else {
  out_dir <- args[1]
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# --- Parameters ---
N <- 100        # individuals
M <- 200        # SNPs
n_pheno <- 3    # phenotypes
chr_counts <- c(80, 70, 50)  # SNPs per chromosome (sum = M)

stopifnot(sum(chr_counts) == M)
n_chr <- length(chr_counts)

# --- 1. Generate SNP data ---
snp_ids <- paste0("rs", seq_len(M))
chromosomes <- rep(seq_len(n_chr), chr_counts)
positions <- unlist(lapply(chr_counts, function(n) sort(sample(1:1e6, n))))
# Genetic distance = 0 (unused for test)
genetic_dist <- rep(0, M)
# Alleles
allele1 <- sample(c("A", "C", "G", "T"), M, replace = TRUE)
allele2 <- sample(c("A", "C", "G", "T"), M, replace = TRUE)
# Ensure they differ
for (i in seq_len(M)) {
  while (allele2[i] == allele1[i]) allele2[i] <- sample(c("A", "C", "G", "T"), 1)
}

# --- 2. Generate genotypes (random 0/1/2 with MAF ~0.3) ---
maf <- runif(M, 0.1, 0.5)
genotypes <- matrix(integer(0), nrow = N, ncol = M)
for (j in seq_len(M)) {
  p <- maf[j]
  # Hardy-Weinberg: P(0)=(1-p)^2, P(1)=2p(1-p), P(2)=p^2
  probs <- c((1 - p)^2, 2 * p * (1 - p), p^2)
  genotypes[, j] <- sample(0:2, N, replace = TRUE, prob = probs)
}

# --- 3. Generate phenotypes (continuous, with some genetic architecture) ---
# Each phenotype: y = sum(beta_j * G_j) + noise
# Choose 10 causal SNPs per phenotype
phenos <- matrix(NA, nrow = N, ncol = n_pheno)
pheno_names <- c("phenoA", "phenoB", "phenoC")

for (k in seq_len(n_pheno)) {
  causal_idx <- sort(sample(seq_len(M), 10))
  causal_effects <- rnorm(10, mean = 0, sd = 0.3)
  # Compute genetic component
  genetic_component <- genotypes[, causal_idx, drop = FALSE] %*% causal_effects
  # Add noise (heritability ~0.4)
  var_g <- var(as.vector(genetic_component))
  noise <- rnorm(N, mean = 0, sd = sqrt(var_g * (1 - 0.4) / 0.4))
  phenos[, k] <- as.vector(genetic_component) + noise
  # Rescale to realistic range (e.g., 3000-5000) without external packages
  r <- range(phenos[, k])
  phenos[, k] <- 3000 + (phenos[, k] - r[1]) / (r[2] - r[1]) * 2000
}

# --- 4. Write FAM file ---
# FID, IID, father, mother, sex, phenotype
fid <- paste0("FAM", sprintf("%03d", seq_len(N)))
iid <- paste0("IND", sprintf("%03d", seq_len(N)))

fam <- data.frame(
  FID = fid,
  IID = iid,
  father = rep(0, N),
  mother = rep(0, N),
  sex = rep(1, N),
  phenotype = rep(-9, N),  # missing pheno in .fam; we'll use multi_pheno_file
  stringsAsFactors = FALSE
)
write.table(fam, file.path(out_dir, "study_sample.fam"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# --- 5. Write BIM file ---
bim <- data.frame(
  chr = chromosomes,
  rsid = snp_ids,
  gdist = genetic_dist,
  pos = positions,
  a1 = allele1,
  a0 = allele2,
  stringsAsFactors = FALSE
)
write.table(bim, file.path(out_dir, "study_sample.bim"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
# Copy as GWAS BIM
file.copy(file.path(out_dir, "study_sample.bim"),
          file.path(out_dir, "gwas.bim"), overwrite = TRUE)

# --- 6. Write BED file (PLINK binary, SNP-major mode) ---
# Format: magic (3 bytes) + SNP-major packed data
# Each SNP: ceil(N/4) bytes, 2-bit genotypes packed 4-per-byte
# Encoding: 00=homoz ref(0), 01=heteroz(1), 10=homoz alt(2), 11=missing

bytes_per_snp <- ceiling(N / 4)
bed_data <- raw(M * bytes_per_snp + 3)

# Magic bytes
bed_data[1] <- as.raw(0x6C)
bed_data[2] <- as.raw(0x1B)
bed_data[3] <- as.raw(0x01)  # SNP-major mode

# Fill genotype data
idx <- 4  # starting byte after magic
for (j in seq_len(M)) {
  for (b in seq_len(bytes_per_snp)) {
    val <- 0L
    for (bb in 0:3) {
      ind <- (b - 1) * 4 + bb + 1
      if (ind <= N) {
        gt <- genotypes[ind, j]
        gt_code <- switch(as.character(gt), "0" = 0L, "1" = 1L, "2" = 2L, 3L)
        val <- val + gt_code * (4L^bb)
      }
    }
    bed_data[idx] <- as.raw(val)
    idx <- idx + 1
  }
}

writeBin(bed_data, file.path(out_dir, "study_sample.bed"))

# --- 7. Generate GWAS summary statistics ---
# For each phenotype, create a "GWAS" by running a simple linear regression
message("Generating summary statistics...")

for (k in seq_len(n_pheno)) {
  y <- phenos[, k]
  ss <- data.frame(
    rsid = snp_ids,
    CHR = chromosomes,
    BP = positions,
    A1 = allele1,
    beta = numeric(M),
    sebeta = numeric(M),
    p = numeric(M),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(M)) {
    fit <- lm(y ~ genotypes[, j])
    coef_summary <- summary(fit)$coef[2, ]
    ss$beta[j] <- coef_summary[1]
    ss$sebeta[j] <- coef_summary[2]
    ss$p[j] <- coef_summary[4]
  }

  out_file <- file.path(out_dir, paste0("gwas_", pheno_names[k], ".txt"))
  write.table(ss, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Wrote: ", out_file)
}

# --- 8. Write combined phenotype file ---
pheno_df <- data.frame(
  FID = fid,
  IID = iid,
  stringsAsFactors = FALSE
)
for (k in seq_len(n_pheno)) {
  pheno_df[[pheno_names[k]]] <- phenos[, k]
}

write.table(pheno_df, file.path(out_dir, "combined_phenotypes.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE)
message("  Wrote: ", file.path(out_dir, "combined_phenotypes.txt"))

# --- 9. Write test config ---
abs_out <- normalizePath(out_dir)
# Assume running from repo root (e.g., Rscript src/generate_test_data.R ...)
detect_repo <- getwd()
# Quick sanity check: if src/generate_test_data.R exists here, it's the repo root
if (!file.exists(file.path(detect_repo, "src/generate_test_data.R"))) {
  detect_repo <- "/path/to/prs_pipeline"
  message("NOTE: Could not auto-detect prs_pipeline root. Set path_repo in test_config.txt")
}

config_content <- sprintf(
'# Test config for multi-phenotype PRS pipeline
# Generated by generate_test_data.R
# Run: sbatch run_single_ancestry_PRS_pipeline.sh -C test_config.txt -c -l -s -P

study_sample="%s/study_sample"
summary_stats_files="%s/gwas_phenoA.txt,%s/gwas_phenoB.txt,%s/gwas_phenoC.txt"
multi_pheno_file="%s/combined_phenotypes.txt"
bim_file_path="%s/gwas.bim"
output_path="%s"
path_repo="%s"
n_total_gwas=100
skip_ss_generation=0
phenotype_info_file="%s/combined_phenotypes.txt"
ncores=4
gwas_pca_eigenvec_file=""

# NOTE: Set path_repo to the absolute path of your prs_pipeline checkout
# path_repo="/home/user/prs_pipeline"

RUN_CT=true
RUN_LDPRED2=true
RUN_LASSOSUM2=true
RUN_PRSice2=true
',
  abs_out,
  abs_out, abs_out, abs_out,
  abs_out,
  abs_out,
  abs_out,
  detect_repo,
  abs_out)

writeLines(config_content, file.path(out_dir, "test_config.txt"))
message("  Wrote: ", file.path(out_dir, "test_config.txt"))

# --- 10. Also generate aligned-format summary stats for test mode ---
# The pipeline's prepare_sumstats.R produces this format.
# A separate config with skip_ss_generation=1 can use these directly.
message("Generating aligned-format summary stats (skip_ss_generation=1)...")
for (k in seq_len(n_pheno)) {
  ss_aligned <- data.frame(
    SNP = snp_ids,
    CHR = chromosomes,
    BP = positions,
    A1 = allele1,
    A2 = allele2,
    beta = as.numeric(NA),
    beta_se = as.numeric(NA),
    P = as.numeric(NA),
    n_eff = rep(N, M),
    stringsAsFactors = FALSE
  )

  y <- phenos[, k]
  for (j in seq_len(M)) {
    fit <- lm(y ~ genotypes[, j])
    cs <- summary(fit)$coef[2, ]
    ss_aligned$beta[j] <- cs[1]
    ss_aligned$beta_se[j] <- cs[2]
    ss_aligned$P[j] <- cs[4]
  }

  out_file <- file.path(out_dir, paste0("aligned_gwas_", pheno_names[k], ".txt"))
  write.table(ss_aligned, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("  Wrote: ", out_file)
}

# --- 11. Write skip-ss config (pre-aligned, skips prepare_sumstats.R) ---
config_skip <- sprintf(
'# Test config for multi-phenotype PRS pipeline (skip_ss_generation=1)
# Generated by generate_test_data.R
# Uses pre-aligned summary stats, skips prepare_sumstats.R

study_sample="%s/study_sample"
summary_stats_files="%s/aligned_gwas_phenoA.txt,%s/aligned_gwas_phenoB.txt,%s/aligned_gwas_phenoC.txt"
multi_pheno_file="%s/combined_phenotypes.txt"
bim_file_path="%s/gwas.bim"
output_path="%s"
path_repo="%s"
n_total_gwas=100
skip_ss_generation=1
phenotype_info_file="%s/combined_phenotypes.txt"
ncores=4
gwas_pca_eigenvec_file=""

RUN_CT=true
RUN_LDPRED2=true
RUN_LASSOSUM2=true
RUN_PRSice2=true
',
  abs_out,
  abs_out, abs_out, abs_out,
  abs_out,
  abs_out,
  abs_out,
  detect_repo,
  abs_out)

writeLines(config_skip, file.path(out_dir, "test_config_skip_ss.txt"))
message("  Wrote: ", file.path(out_dir, "test_config_skip_ss.txt"))

message("\nDone! Test data generated in: ", abs_out)
message("\nConfigs:")
message("  (1) Full pipeline (runs prepare_sumstats.R):")
message("       sbatch run_single_ancestry_PRS_pipeline.sh -C ", abs_out, "/test_config.txt -c -l -s -P")
message("  (2) Skip SS gen (pre-aligned, for quick testing):")
message("       sbatch run_single_ancestry_PRS_pipeline.sh -C ", abs_out, "/test_config_skip_ss.txt -c -l -s -P")
message("\nOr via singularity:")
message("  sbatch sandbox_singularity_runner.sh -C ", abs_out, "/test_config_skip_ss.txt -c -l -s -P")
