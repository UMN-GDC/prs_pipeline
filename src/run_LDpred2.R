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
library(data.table)
options(bigstatsr.check.parallel = FALSE)
options(bigstatsr.check.args = FALSE)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='LDpred2 Pipeline for Polygenic Risk Scores')

parser$add_argument("--anc_bed", type="character", help="Path to the plink bed file")
parser$add_argument("--rds", type="character", help="Path to the .rds file")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to the .bim file")

# Both are now optional; logic below handles various scenarios
parser$add_argument("--ss", type="character", required=TRUE, help="Path to summary statistics")
parser$add_argument("--beta_se", type="character", help="Optional: separate table with beta_se/SE")
parser$add_argument("--afreq", type="character", help="Path to PLINK2 .afreq file (bypasses snp_MAF on genotype matrix)")
parser$add_argument("--pheno", type="character", help="Path to phenotype file (FID, IID, phenotype). Overrides .fam affection column.")

parser$add_argument("--h2", type="numeric", default=0.4, help="Assumed heritability")
parser$add_argument("--n_val", type="integer", default=49, help="Number of samples for validation")
parser$add_argument("--out", type="character", default="ldpred2_out", help="Prefix for output files")
parser$add_argument("--ncores", type="integer", default=1, help="Number of CPU cores (default: 1)")
parser$add_argument("--ld-cache-dir", type="character", help="Directory to cache/reuse per-chromosome LD matrices")
parser$add_argument("--ld-matrix-dir", type="character", help="Directory with pre-computed LD matrix (from generate_ld_matrix.R)")

args <- parser$parse_args()

# Wrapper to capture full state on any error
tryCatch({

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
CHR_id <- obj.bigSNP$map$chromosome
POS_id <- obj.bigSNP$map$physical.pos
y      <- obj.bigSNP$fam$affection
NCORES <- args$ncores

# Override phenotype from external file if provided
if (!is.null(args$pheno)) {
  message("Reading phenotype from: ", args$pheno)
  pheno <- fread2(args$pheno, header = TRUE)
  y <- pheno[[3]]
  if (length(y) != length(obj.bigSNP$fam$affection)) {
    stop("Phenotype file has ", length(y), " rows, but the .rds file has ",
         length(obj.bigSNP$fam$affection), " individuals.")
  }
}

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

# Merge with BIM to get a0 (reference allele)
sumstats <- sumstats %>%
  inner_join(bim.file, by = c("chr", "pos")) %>%
  mutate(a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1)) %>%
  select(chr, pos, a1, a0, beta, beta_se, n_eff)
message("SNPs matching between Summary Stats and BIM: ", nrow(sumstats))

# --- 4. PREPARE FOR LDPRED2 ---

if (!is.null(args$ld_matrix_dir)) {
  message("Aligning sumstats with pre-computed LD matrix map...")
  map <- readRDS(file.path(args$ld_matrix_dir, "map.rds"))
} else {
  map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
}
df_beta <- snp_match(sumstats, map, join_by_pos = TRUE)

# When using pre-computed LD matrix, remap _NUM_ID_ from map index to G column index
if (!is.null(args$ld_matrix_dir)) {
  g_idx <- readRDS(file.path(args$ld_matrix_dir, "g_idx.rds"))
  df_beta$`_LOCAL_ID_` <- df_beta$`_NUM_ID_`
  df_beta$`_NUM_ID_` <- g_idx[df_beta$`_NUM_ID_`]
}

POS2 <- obj.bigSNP$map$genetic.dist
ind.row <- rows_along(G)

if (!is.null(args$afreq)) {
  message("Reading allele frequencies from: ", args$afreq)
  afreq <- fread2(args$afreq)
  m <- match(df_beta$rsid, afreq$ID)
  maf <- pmin(afreq$ALT_FREQS[m], 1 - afreq$ALT_FREQS[m])
  maf[is.na(m)] <- NA
} else {
  message("Computing MAF from genotype matrix...")
  maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = NCORES)
}

maf_thr <- 1 / sqrt(length(ind.row))
message("--- MAF Filtering Diagnostics ---")
message("   Total variants matched: ", length(maf))
message("   Variants with NA MAF (missing data): ", sum(is.na(maf)))
message("   Variants failing threshold (<= ", round(maf_thr, 4), "): ", sum(maf <= maf_thr, na.rm=TRUE))
message("   Variants passing threshold: ", sum(maf > maf_thr, na.rm=TRUE))

df_beta <- df_beta[maf > maf_thr & !is.na(maf), ]

# Restrict to autosomes 1:22
CHRS <- 1:22
df_beta <- df_beta[df_beta$chr %in% CHRS, ]
message("   Autosomes 1-22 after MAF filter: ", nrow(df_beta))

if (nrow(df_beta) == 0) {
  stop("No variants remain after MAF + autosome filtering.")
}

# --- 5. LD MATRIX ---

if (!is.null(args$ld_cache_dir)) {
  dir.create(args$ld_cache_dir, showWarnings = FALSE, recursive = TRUE)
}

tmp <- tempfile(tmpdir = "temp_ld")
dir.create("temp_ld", showWarnings = FALSE)
corr <- NULL
ld <- NULL
keep_idx <- logical(nrow(df_beta))

if (!is.null(args$ld_matrix_dir)) {
  # Load pre-computed LD matrix (sumstats-independent)
  message("Loading pre-computed LD matrix from: ", args$ld_matrix_dir)
  ld_map <- readRDS(file.path(args$ld_matrix_dir, "map.rds"))
  message("  ld_map rows: ", nrow(ld_map), ", cols: ", paste(colnames(ld_map), collapse=","))

  for (chr in CHRS) {
    ind.chr <- which(df_beta$chr == chr)
    if (length(ind.chr) < 2) next

    map_chr_idx <- which(ld_map$chr == chr)
    message("  Chr", chr, ": df_beta has ", length(ind.chr), " SNPs, ld_map has ", length(map_chr_idx), " SNPs")
    local_idx <- match(df_beta$`_LOCAL_ID_`[ind.chr], map_chr_idx)
    message("    _LOCAL_ID_ range: [", min(df_beta$`_LOCAL_ID_`[ind.chr], na.rm=TRUE), ", ", max(df_beta$`_LOCAL_ID_`[ind.chr], na.rm=TRUE), "], map_chr_idx range: [", min(map_chr_idx), ", ", max(map_chr_idx), "]")
    bad <- which(is.na(local_idx))
    if (length(bad) > 0) {
      local_idx <- local_idx[-bad]
      ind.chr <- ind.chr[-bad]
    }
    if (length(local_idx) < 2) next

    keep_idx[ind.chr] <- TRUE
    corr0_full <- readRDS(file.path(args$ld_matrix_dir, sprintf("chr%d_corr.rds", chr)))
    corr0 <- corr0_full[local_idx, local_idx]

    if (any(is.na(corr0))) {
      message("  Chr", chr, ": zeroing ", sum(is.na(corr0)), " NA correlations")
      corr0[is.na(corr0)] <- 0
    }
    # Add small ridge to ensure positive definiteness for the Gibbs sampler
    corr0 <- corr0 + diag(ncol(corr0)) * 1e-5

    if (is.null(corr)) {
      ld <- Matrix::colSums(corr0^2)
      corr <- as_SFBM(corr0, tmp, compact = TRUE)
    } else {
      ld <- c(ld, Matrix::colSums(corr0^2))
      corr$add_columns(corr0, nrow(corr))
    }
  }
} else {
  # Check if per-chromosome cache exists
  use_cache <- !is.null(args$ld_cache_dir) &&
    all(vapply(CHRS, function(chr) {
      file.exists(file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
    }, logical(1)))

  if (use_cache) {
    message("Loading cached LD matrices from: ", args$ld_cache_dir)
    for (chr in CHRS) {
      ind.chr <- which(df_beta$chr == chr)
      if (length(ind.chr) < 2) next
      corr0 <- readRDS(file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
      if (any(is.na(corr0))) {
        message("  Chr", chr, ": zeroing ", sum(is.na(corr0)), " NA correlations (cached)")
        corr0[is.na(corr0)] <- 0
      }
      corr0 <- corr0 + diag(ncol(corr0)) * 1e-5
      keep_idx[ind.chr] <- TRUE
      if (is.null(corr)) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
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
      keep_idx[ind.chr] <- TRUE

      corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3/1000, infos.pos = POS2[ind.chr2], ncores = NCORES)

      if (any(is.na(corr0))) {
        na_count <- sum(is.na(corr0))
        message("  Chr", chr, ": repairing ", na_count, " NA correlations")
        corr0[is.na(corr0)] <- 0
      }
      corr0 <- corr0 + diag(ncol(corr0)) * 1e-5

      if (!is.null(args$ld_cache_dir)) {
        saveRDS(corr0, file.path(args$ld_cache_dir, paste0("chr", chr, "_corr.rds")))
        message("  Cached chr", chr, " (", length(ind.chr2), " SNPs)")
      }

      if (is.null(corr)) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
    }
  }
}

# Align df_beta with corr
df_beta <- df_beta[keep_idx, ]
message("Variants in LD matrix: ", nrow(df_beta))

if (nrow(df_beta) != ncol(corr)) {
  stop("Mismatch between df_beta rows (", nrow(df_beta),
       ") and corr cols (", ncol(corr), ") after LD construction.")
}

# --- 6. RUN MODELS ---

# Split validation/test sets
n_total <- nrow(G)
n_val <- min(args$n_val, floor(n_total / 3))
if (n_val < 2) stop("Too few samples (", n_total, ") for validation split.")
ind.val  <- sample(n_total, n_val)
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

# Diagnostic: report how many grid points have NA betas
na_beta_cols <- which(apply(beta_grid, 2, function(col) all(is.na(col))))
non_na_beta_cols <- which(apply(beta_grid, 2, function(col) any(!is.na(col))))
message("Beta_grid: ", length(non_na_beta_cols), "/", ncol(beta_grid), " columns have non-NA estimates")
if (length(na_beta_cols) > 0) {
  message("  NA columns (", length(na_beta_cols), "): first few params — ",
    paste(capture.output(print(params[head(na_beta_cols, 3), ])), collapse = " | "))
}

pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])

# Tuning on Validation set
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  if (sd(x, na.rm = TRUE) == 0) return(NA)
  summary(lm(y[ind.val] ~ x))$coef[2, 3] # Use t-stat as score
})

# Save Grid Plot
p <- ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  geom_point() + geom_line() + scale_x_log10() +
  facet_wrap(~ sparse) + labs(title="LDpred2 Grid Search Tuning")
ggsave(paste0(args$out, "_grid_plot.png"), p)

# Best Grid Prediction
best_grid_idx <- which.max(params$score)
if (length(best_grid_idx) == 0 || is.na(best_grid_idx) || all(is.na(params$score))) {
  message("WARNING: All grid model scores are NA. Skipping grid PRS.")
  pred_grid_best <- rep(NA_real_, length(ind.test))
  r2_grid <- NA_real_
  prs_grid_all <- rep(NA_real_, n_total)
} else {
  pred_grid_best <- big_prodVec(G, beta_grid[, best_grid_idx], ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
  r2_grid <- pcor(pred_grid_best, y[ind.test], NULL)

  # Calculate scores for EVERYONE (not just the test set)
  prs_grid_all <- pred_grid[, best_grid_idx]
}

# --- 7. OUTPUT RESULTS ---

# Recalculate Infinitesimal for everyone (removing the ind.row restriction)
message("G is :", dim(G), " and the G matrix head (10x10):", G[1:10, 1:10])
message("Head on beta inf:", head(beta_inf))

prs_inf_all <- big_prodVec(G, beta_inf, ind.col = df_beta[["_NUM_ID_"]])

# 2. Create the individual-level table
# obj.bigSNP$fam contains family.ID (FID) and sample.ID (IID)
prs_report <- data.frame(
  FID      = obj.bigSNP$fam$family.ID,
  IID      = obj.bigSNP$fam$sample.ID,
  PRS_inf  = prs_inf_all,
  PRS_grid = prs_grid_all
)

# 3. Save the individual scores
fwrite(prs_report, paste0(args$out, "_individual_scores.txt"), sep = "\t", row.names = FALSE)

# 4. Save the performance summary
results <- data.frame(
  Method = c("Infinitesimal", "Grid"),
  R2 = c(r2_inf^2, r2_grid^2)
)
write.csv(results, paste0(args$out, "_performance.csv"), row.names = FALSE)

message("Success! Individual scores saved to: ", paste0(args$out, "_individual_scores.txt"))

}, error = function(e) {
  message("\n========== ERROR IN run_LDpred2.R ==========")
  message("Condition: ", conditionMessage(e))
  message("Call stack:")
  for (i in seq_len(sys.nframe())) {
    call <- sys.call(i)
    if (!is.null(call)) message("  ", i, ": ", deparse(call)[1])
  }
  dump.frames("ldpred2_dump", to.file = TRUE)
  message("Dumped to ldpred2_dump.rda")
  stop(e)
})
