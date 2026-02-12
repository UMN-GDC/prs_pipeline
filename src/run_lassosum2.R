#!/usr/bin/env Rscript

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='lassosum2 Pipeline for Polygenic Risk Scores')
parser$add_argument("--anc_bed", type="character", help="Path to the plink .bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--ss", type="character", required=TRUE, help="Path to GWAS summary statistics")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")
parser$add_argument("--beta_se", type="character", required=TRUE, help="Path to table with rsid, chr, and beta_se")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation (default 49)")
parser$add_argument("--seed", type="integer", default=1, help="Seed for validation/test split (default 1)")
parser$add_argument("--out", type="character", default="lassosum_out", help="Prefix for output files")

args <- parser$parse_args()

# 2. DATA LOADING LOGIC
if (!is.null(args$anc_bed)) {
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)
  if (!file.exists(rds_path)) {
    message("Converting .bed to .rds...")
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path
}

if (is.null(args$rds) || !file.exists(args$rds)) {
  stop("Error: Valid .rds file path required. Provide --anc_bed or a correct --rds path.")
}

obj.bigSNP <- snp_attach(args$rds)
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()

# 3. SUMMARY STATS PREPARATION
qassoc <- read.table(args$ss, header = TRUE, sep = "", stringsAsFactors = FALSE) %>% 
  filter(TEST == "ADD") %>% 
  rename(rsid = SNP, pos = BP, chr = CHR)

bim.file <- bigreadr::fread2(args$bim, select = c(1, 4, 5, 6))
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")

beta.se.tab <- read.table(args$beta_se, header = TRUE, sep = "", stringsAsFactors = FALSE) %>% 
  select(rsid, chr, beta_se)

sumstats <- qassoc %>% 
  transmute(rsid, chr, pos, a1=A1, beta=BETA, n_eff=NMISS) %>% 
  inner_join(bim.file, by=c("chr", "pos")) %>%
  mutate(a0=ifelse(a1==bim.a1, bim.a0, bim.a1)) %>%
  inner_join(beta.se.tab, by=c("chr", "rsid")) %>%
  select(chr, pos, a1, a0, beta, beta_se, n_eff)

map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

# MAF Filtering
maf <- snp_MAF(G, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
df_beta <- df_beta[maf > (1 / sqrt(nrow(G))), ]

# 4. LD MATRIX COMPUTATION
tmp <- tempfile(tmpdir = "temp_ld_lassosum")
dir.create("temp_ld_lassosum", showWarnings = FALSE)
POS2 <- obj.bigSNP$map$genetic.dist
corr <- NULL

for (chr in 1:22) {
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  if (length(ind.chr2) < 2) next
    
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)
    
  if (is.null(corr)) {
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr0, nrow(corr))
  }
}

# 5. LASSOSUM2 MODELING
message("Running lassosum2 grid search...")
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")

# Split validation/test
set.seed(args$seed)
ind.val  <- sample(nrow(G), args$n_val)
ind.test <- setdiff(rows_along(G), ind.val)

# Predict on validation set for tuning
pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]])
params2$score <- apply(pred_grid2[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  # Using the Z-score (column 3 of coefficients)
  summary(lm(y[ind.val] ~ x))$coef["x", 3]
})

# 6. RESULTS & VISUALIZATION
# Plot tuning results
p <- ggplot(params2, aes(x = lambda, y = score, color = as.factor(delta))) +
  theme_bigstatsr() +
  geom_point() + geom_line() +
  scale_x_log10() +
  labs(title = "lassosum2 Tuning", y = "Z-Score", color = "delta")
ggsave(paste0(args$out, "_lassosum_plot.png"), p)

# Best Model Performance
best_idx <- which.max(params2$score)
best_beta <- beta_lassosum2[, best_idx]
pred_test <- big_prodVec(G, best_beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
test_r2 <- pcor(pred_test, y[ind.test], NULL)^2

# Save performance and best parameters
write.csv(params2, paste0(args$out, "_grid_params.csv"), row.names = FALSE)
cat(paste("Best Test R2:", test_r2, "\n"), file = paste0(args$out, "_final_res.txt"))

# Clean up temp files
unlink("temp_ld_lassosum", recursive = TRUE)
message("Success! Best Test R2: ", test_r2)