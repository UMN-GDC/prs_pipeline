#!/usr/bin/env Rscript

# Control BLAS threads to avoid nested parallelism with bigstatsr's multicore
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(bigreadr)

options(bigstatsr.check.parallel = FALSE)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='lassosum2 Pipeline for Polygenic Risk Scores')
parser$add_argument("--anc_bed", type="character", help="Path to the plink .bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--ss", type="character", required=TRUE, help="Path to GWAS summary statistics")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")
parser$add_argument("--afreq", type="character", help="Path to PLINK2 .afreq file (bypasses snp_MAF on genotype matrix)")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--seed", type="integer", default=1, help="Seed for validation/test split")
parser$add_argument("--out", type="character", default="lassosum_out", help="Prefix for output files")
parser$add_argument("--ncores", type="integer", default=1, help="Number of CPU cores (default: 1)")
parser$add_argument("--ld-cache-dir", type="character", help="Directory to cache/reuse per-chromosome LD matrices")

args <- parser$parse_args()

# --- 2. DATA LOADING ---
if (!is.null(args$anc_bed)) {
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
  if (!file.exists(rds_path)) {
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path
}

obj.bigSNP <- snp_attach(args$rds)
G      <- obj.bigSNP$genotypes
y      <- obj.bigSNP$fam$affection
NCORES <- args$ncores

# --- 3. SUMMARY STATS PREPARATION (Adjusted for your new columns) ---
# Use fread2 to handle the 9-column file
qassoc <- bigreadr::fread2(args$ss) 

# Apply the same mapping logic from the LDpred2 script
sumstats_mapped <- qassoc %>%
  rename(
    rsid  = any_of(c("SNP", "rsid", "ID")),
    pos   = any_of(c("BP", "pos", "position")),
    chr   = any_of(c("CHR", "chr", "chromosome")),
    a1    = any_of(c("A1", "a1", "allele1")),
    beta  = any_of(c("BETA", "beta", "eff")),
    n_eff = any_of(c("n_eff", "NMISS", "N"))
  )

if ("TEST" %in% colnames(sumstats_mapped)) {
  sumstats_mapped <- sumstats_mapped %>% filter(TEST == "ADD")
}

bim.file <- bigreadr::fread2(args$bim, select = c(1, 4, 5, 6))
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")

# Handle SE logic exactly as in the SS processing
if (!is.null(args$beta_se)) {
  beta.se.tab <- bigreadr::fread2(args$beta_se) %>% 
    rename(rsid = any_of(c("SNP", "rsid")), 
           beta_se = any_of(c("SE", "beta_se"))) %>%
    select(rsid, beta_se)
  
  sumstats <- sumstats_mapped %>% inner_join(beta.se.tab, by = "rsid")
} else {
  sumstats <- sumstats_mapped %>% rename(beta_se = any_of(c("SE", "beta_se")))
}

# Final alignment
sumstats <- sumstats %>%
  inner_join(bim.file, by = c("chr", "pos")) %>%
  mutate(a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1)) %>%
  select(chr, pos, a1, a0, beta, beta_se, n_eff)

map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

# MAF Filtering
if (!is.null(args$afreq)) {
  message("Reading allele frequencies from: ", args$afreq)
  afreq <- fread2(args$afreq)
  m <- match(df_beta$rsid, afreq$ID)
  maf <- pmin(afreq$ALT_FREQS[m], 1 - afreq$ALT_FREQS[m])
  maf[is.na(m)] <- NA
} else {
  message("Computing MAF from genotype matrix...")
  maf <- snp_MAF(G, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
}
df_beta <- df_beta[maf > (1 / sqrt(nrow(G))) & !is.na(maf), ]

# --- 4. LD MATRIX COMPUTATION ---
if (!is.null(args$ld_cache_dir)) {
  dir.create(args$ld_cache_dir, showWarnings = FALSE, recursive = TRUE)
}

tmp <- tempfile(tmpdir = "temp_ld_lassosum")
dir.create("temp_ld_lassosum", showWarnings = FALSE)
POS2 <- obj.bigSNP$map$genetic.dist
corr <- NULL

CHRS <- 1:22

# Check if all chromosomes are cached
use_cache <- !is.null(args$ld_cache_dir) &&
  all(vapply(CHRS, function(chr) {
    file.exists(file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
  }, logical(1)))

if (use_cache) {
  message("Loading pre-computed LD matrices from: ", args$ld_cache_dir)
  for (chr in CHRS) {
    ind.chr <- which(df_beta$chr == chr)
    if (length(ind.chr) < 2) next
    corr0 <- readRDS(file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
    if (is.null(corr)) {
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr0, nrow(corr))
    }
  }
} else {
  message("Computing LD matrices from genotype data...")
  for (chr in CHRS) {
    ind.chr <- which(df_beta$chr == chr)
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    if (length(ind.chr2) < 2) next

    message(paste("Processing Chromosome:", chr, "| SNPs:", length(ind.chr2)))

    corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)

    if (!is.null(args$ld_cache_dir)) {
      saveRDS(corr0, file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
      message("  Cached chr", chr, " (", length(ind.chr2), " SNPs)")
    }

    if (is.null(corr)) {
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr0, nrow(corr))
    }
  }
}

# --- 5. LASSOSUM2 MODELING ---
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")

set.seed(args$seed)
ind.val  <- sample(nrow(G), args$n_val)
ind.test <- setdiff(rows_along(G), ind.val)

pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]])
params2$score <- apply(pred_grid2[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef[2, 3]
})

# --- 6. RESULTS ---
p <- ggplot(params2, aes(x = lambda, y = score, color = as.factor(delta))) +
  theme_bigstatsr() + geom_point() + geom_line() + scale_x_log10() +
  labs(title = "lassosum2 Tuning", y = "Z-Score", color = "delta")
ggsave(paste0(args$out, "_lassosum_plot.png"), p)

best_idx <- which.max(params2$score)
best_beta <- beta_lassosum2[, best_idx]
pred_test <- big_prodVec(G, best_beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
test_r2 <- pcor(pred_test, y[ind.test], NULL)^2

write.csv(params2, paste0(args$out, "_grid_params.csv"), row.names = FALSE)
cat(paste("Best Test R2:", test_r2, "\n"), file = paste0(args$out, "_final_res.txt"))

unlink("temp_ld_lassosum", recursive = TRUE)
message(paste0("Success! Best Test R2: ", round(test_r2, 4)))

# --- 7. Save RESULTS ---
# pred_grid2 contains rows = individuals, columns = grid points
full_results <- data.frame(
  FID = obj.bigSNP$fam[, 1],
  IID = obj.bigSNP$fam[, 2],
  pred_grid2
)

# Rename columns to match the row index of the grid_params.csv for easy lookup
colnames(full_results)[3:ncol(full_results)] <- paste0("grid_", 1:nrow(params2))

# Write to a compressed CSV or text file (it can be large)
bigreadr::fwrite2(full_results, paste0(args$out, "_full_predictions.csv"))

# --- 8. Simplify the process of getting SCORE column ---
# Find the index of the best performing model
best_idx <- which.max(params2$score)
best_col_name <- paste0("grid_", best_idx)

# Extract only FID, IID, and the Best PRS
final_prs_df <- full_results[, c("FID", "IID", best_col_name)]

# Rename the column for clarity
colnames(final_prs_df)[3] <- "Best_PRS_Score"

# Save the finalized single-score file
bigreadr::fwrite2(final_prs_df, paste0(args$out, "_final_best_prs.csv"))

message(paste("Selected", best_col_name, "as the best model based on validation score."))