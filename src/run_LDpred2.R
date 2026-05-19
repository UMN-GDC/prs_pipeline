#!/usr/bin/env Rscript

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)
library(bigreadr)
library(data.table)
options(bigstatsr.check.parallel = FALSE)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='LDpred2 Pipeline for Polygenic Risk Scores')

parser$add_argument("--anc_bed", type="character", help="Path to the plink bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")

# Both are now optional; logic below handles various scenarios
parser$add_argument("--ss", type="character", required=TRUE, help="Path to summary statistics")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")

parser$add_argument("--h2", type="numeric", default=0.4, help="Assumed heritability")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--out", type="character", default="ldpred2_out", help="Prefix for output files")

args <- parser$parse_args()

# --- 2. DATA LOADING & CONVERSION ---

if (!is.null(args$anc_bed)) {
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
  if (!file.exists(rds_path)) {
    message("Converting .bed to .rds format...")
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path 
}

if (is.null(args$rds)) stop("Error: You must provide either --anc_bed or --rds.")

obj.bigSNP <- snp_attach(args$rds)
G      <- obj.bigSNP$genotypes

NCORES <- 1
# NEW: Check for missing values and impute them on the fly
if (any(is.na(G[]))) {
  message("Warning: Target genotype matrix contains missing values (NAs).")
  message("Imputing missing calls using the most frequent allele (mode)...")
  G <- snp_fastImputeSimple(G, method = "mode", ncores = NCORES)
}

CHR_id <- obj.bigSNP$map$chromosome
POS_id <- obj.bigSNP$map$physical.pos
y      <- obj.bigSNP$fam$affection


# Load BIM for alignment
bim.file <- fread2(args$bim, select = c(1, 4, 5, 6))
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")

# --- 3. SUMMARY STATISTICS PROCESSING ---

message("Reading summary statistics from: ", args$ss)
ss_raw <- fread2(args$ss)

# Map common GWAS headers to internal names
sumstats <- ss_raw %>%
  rename(
    rsid  = any_of(c("SNP", "rsid", "ID")),
    pos   = any_of(c("BP", "pos", "position")),
    chr   = any_of(c("CHR", "chr", "chromosome")),
    a1    = any_of(c("A1", "a1", "allele1")),
    beta  = any_of(c("BETA", "beta", "eff")),
    n_eff = any_of(c("n_eff", "NMISS", "N"))
  )

# Only filter by TEST if the column actually exists (e.g., PLINK outputs)
if ("TEST" %in% colnames(sumstats)) {
  sumstats <- sumstats %>% filter(TEST == "ADD")
}

# Logic for SE (Standard Error)
if (!is.null(args$beta_se)) {
  message("Merging with external SE file: ", args$beta_se)
  beta_se_tab <- fread2(args$beta_se) %>% 
    rename(rsid = any_of(c("SNP", "rsid")), 
           beta_se = any_of(c("SE", "beta_se", "beta_se_tab"))) %>%
    select(rsid, beta_se)
  
  sumstats <- sumstats %>% inner_join(beta_se_tab, by = "rsid")
} else {
  message("Extracting SE from primary summary stat file...")
  sumstats <- sumstats %>% rename(beta_se = any_of(c("SE", "beta_se")))
}

# Ensure required columns exist after mapping
req_cols <- c("chr", "pos", "a1", "beta", "beta_se", "n_eff")
missing <- setdiff(req_cols, colnames(sumstats))
if (length(missing) > 0) {
  stop("Missing required columns in summary stats: ", paste(missing, collapse = ", "))
}

# Keep rsid in the select statement!
sumstats <- sumstats %>%
  inner_join(bim.file, by = c("chr", "pos")) %>%
  mutate(a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1)) %>%
  select(chr, rsid, pos, a1, a0, beta, beta_se, n_eff)

# DIAGNOSTIC: Check if the data actually matched
message("SNPs matching between Summary Stats and BIM: ", nrow(sumstats))
if (nrow(sumstats) == 0) {
  stop("CRITICAL ERROR: 0 SNPs matched! Check your chromosome column formatting (e.g., '01' vs '1').")
}

message("SNPs matching between Summary Stats and BIM: ", nrow(sumstats))

# --- 4. PREPARE FOR LDPRED2 ---
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

# Fallback window selection: Use physical position if genetic distance map is missing (all zeros)
if (all(obj.bigSNP$map$genetic.dist == 0)) {
  message("Genetic distance is all zeros. Falling back to physical positions (kb) with 500kb window.")
  POS2 <- obj.bigSNP$map$physical.pos / 1000
  ld_window_size <- 500
} else {
  POS2 <- obj.bigSNP$map$genetic.dist
  ld_window_size <- 3/1000
}

ind.row <- rows_along(G)
maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 0.0001 # Adjust this to smaller numbers to be more permissive

# Diagnostic metrics
message("--- MAF Filtering Diagnostics ---")
message("  Total variants matched: ", length(maf))
message("  Variants with NA MAF (missing data): ", sum(is.na(maf)))
message("  Variants failing threshold (<= ", round(maf_thr, 4), "): ", sum(maf <= maf_thr, na.rm=TRUE))
message("  Variants passing threshold: ", sum(maf > maf_thr, na.rm=TRUE))

# CRITICAL FIX 1: Wrap in which() to prevent generating trailing NA rows
df_beta <- df_beta[which(maf > maf_thr), ]

if (nrow(df_beta) == 0) {
  stop("CRITICAL ERROR: 0 variants passed the MAF threshold filter.")
}

# --- 5. COMPUTE LD MATRIX ---

tmp <- tempfile(tmpdir = "temp_ld")
dir.create("temp_ld", showWarnings = FALSE)
corr <- NULL
ld <- NULL

kept_df_indices <- c()

for (chr in 1:22) {
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  if (length(ind.chr2) < 2) {
    if (length(ind.chr2) == 1) {
      message(paste("Skipping Chromosome:", chr, "| Not enough SNPs (1)"))
    }
    next
  }
  
  message(paste("Processing Chromosome:", chr, "| SNPs:", length(ind.chr2)))
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = ld_window_size, infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (is.null(corr)) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
  
  kept_df_indices <- c(kept_df_indices, ind.chr)
}

# CRITICAL FIX 3: Subset df_beta down to match the exact dimensions of the compiled corr matrix
if (length(kept_df_indices) == 0) {
  stop("CRITICAL ERROR: No chromosomes contained enough valid SNPs to build an LD Matrix.")
}
df_beta <- df_beta[kept_df_indices, ]

# --- 6. RUN MODELS ---

# Split validation/test sets
ind.val  <- sample(nrow(G), args$n_val)
ind.test <- setdiff(rows_along(G), ind.val)

# Infinitesimal Model
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = args$h2)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta$`_NUM_ID_`)
r2_inf   <- pcor(pred_inf, y[ind.test], NULL)

# Grid Model
h2_seq <- round(args$h2 * c(0.7, 1, 1.4), 4)
p_seq  <- signif(seq_log(1e-4, 1, length.out = 10), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

# Tuning on Validation set
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef[2, 3] # Use t-stat as score
})

# Save Grid Plot
p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  geom_point() + geom_line() + scale_x_log10() +
  facet_wrap(~ sparse) + labs(title="LDpred2 Grid Search Tuning")
ggsave(paste0(args$out, "_grid_plot.png"), p)

# Best Grid Prediction
best_grid_idx <- which.max(params$score)
pred_grid_best <- big_prodVec(G, beta_grid[, best_grid_idx], ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
r2_grid <- pcor(pred_grid_best, y[ind.test], NULL)

# --- 7. OUTPUT RESULTS ---

# 1. Calculate scores for EVERYONE 
prs_grid_all <- pred_grid[, best_grid_idx]

# Recalculate Infinitesimal for everyone 
prs_inf_all <- big_prodVec(G, beta_inf, ind.col = df_beta[["_NUM_ID_"]])

# 2. Create individual table
prs_report <- data.frame(
  FID      = obj.bigSNP$fam$family.ID,
  IID      = obj.bigSNP$fam$sample.ID,
  PRS_inf  = prs_inf_all,
  PRS_grid = prs_grid_all
)

# 3. Save individual scores
fwrite(prs_report, paste0(args$out, "_individual_scores.txt"), sep = "\t", row.names = FALSE)

# 4. Save performance summary
results <- data.frame(
  Method = c("Infinitesimal", "Grid"),
  R2 = c(r2_inf^2, r2_grid^2)
)
write.csv(results, paste0(args$out, "_performance.csv"), row.names = FALSE)

message("Success! Individual scores saved to: ", paste0(args$out, "_individual_scores.txt"))