#!/usr/bin/env Rscript

library(argparse)
library(bigsnpr)
library(ggplot2)
library(dplyr)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='LDpred2 Pipeline for Polygenic Risk Scores')
# Changed these to required=FALSE (default)
parser$add_argument("--anc_bed", type="character", help="Path to the plink bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--ss", type="character", required=TRUE, help="Path to summary statistics (assoc.linear)")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")
parser$add_argument("--beta_se", type="character", required=TRUE, help="Path to table with beta_se")
parser$add_argument("--h2", type="numeric", default=0.4, help="Assumed heritability")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--out", type="character", default="ldpred2_out", help="Prefix for output files")

args <- parser$parse_args()

# --- DATA LOADING LOGIC ---

# 1. If anc_bed is provided, we need to make sure an RDS exists
if (!is.null(args$anc_bed)) {
  # Determine the RDS name: use provided --rds or default to the bed filename
  rds_path <- if (!is.null(args$rds)) args$rds else sub("\\.bed$", ".rds", args$anc_bed)

  # Only run snp_readBed if the RDS file doesn't already exist
  if (!file.exists(rds_path)) {
    message("Converting .bed to .rds format...")
    snp_readBed(args$anc_bed)
  }
  args$rds <- rds_path # Ensure args$rds is set for the snp_attach step
}

# 2. Check if we have an RDS path by now
if (is.null(args$rds)) {
  stop("Error: You must provide either --anc_bed (to convert) or --rds (to attach).")
}

# 3. Attach the bigSNP object
obj.bigSNP <- snp_attach(args$rds)
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()

# Read and clean summary stats
qassoc <- read.table(args$ss, header = TRUE, sep = "", stringsAsFactors = FALSE) %>% 
  filter(TEST == "ADD") %>% 
  rename(rsid = SNP, pos = BP, chr = CHR)

bim.file <- bigreadr::fread2(args$bim, select = c(1, 4, 5, 6))
colnames(bim.file) <- c("chr", "pos", "bim.a1", "bim.a0")

beta.se.tab <- read.table(args$beta_se, header = TRUE, sep = "", stringsAsFactors = FALSE) %>% 
  select(rsid, chr, beta_se)

# Merge everything
sumstats <- qassoc %>% 
  transmute(rsid, chr, pos, a1=A1, beta=BETA, n_eff=NMISS) %>% 
  inner_join(bim.file, by=c("chr", "pos")) %>%
  mutate(a0=ifelse(a1==bim.a1, bim.a0, bim.a1)) %>%
  inner_join(beta.se.tab, by=c("chr", "rsid")) %>%
  select(chr, pos, a1, a0, beta, beta_se, n_eff)

# 3. PREPARE FOR LDPRED2
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

POS2 <- obj.bigSNP$map$genetic.dist
ind.row <- rows_along(G)
maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
maf_thr <- 1 / sqrt(length(ind.row))
df_beta <- df_beta[maf > maf_thr, ]

# 4. COMPUTE LD MATRIX
tmp <- tempfile(tmpdir = "temp_ld")
dir.create("temp_ld", showWarnings = FALSE)
corr <- NULL
ld <- NULL

for (chr in 1:22) {
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  if (length(ind.chr2) < 2) next

  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)

  if (is.null(corr)) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

# 5. RUN MODELS
ind.val  <- sample(nrow(G), args$n_val)
ind.test <- setdiff(rows_along(G), ind.val)

# Infinitesimal Model
beta_inf <- snp_ldpred2_inf(as_SFBM(snp_cor(G, ind.col = df_beta$`_NUM_ID_`, size = 3/1000, ncores=NCORES)), 
                            df_beta, h2 = args$h2)
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
  summary(lm(y[ind.val] ~ x))$coef[2, 3]
})

# Save Grid Plot
p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  geom_point() + geom_line() + scale_x_log10() +
  facet_wrap(~ sparse) + labs(title="LDpred2 Grid Search")
ggsave(paste0(args$out, "_grid_plot.png"), p)

# Best Grid Prediction
best_grid_idx <- which.max(params$score)
pred_grid_best <- big_prodVec(G, beta_grid[, best_grid_idx], ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
r2_grid <- pcor(pred_grid_best, y[ind.test], NULL)

# 6. OUTPUT RESULTS
results <- data.frame(
  Method = c("Infinitesimal", "Grid"),
  R2 = c(r2_inf^2, r2_grid^2)
)
write.csv(results, paste0(args$out, "_performance.csv"), row.names = FALSE)

message("Success! Results saved with prefix: ", args$out)