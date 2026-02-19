#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(bigreadr)

# 1. SET UP ARGUMENT PARSER
parser <- ArgumentParser(description='Clean GWAS summary stats and calculate n_eff')
parser$add_argument("--input", type="character", required=TRUE, help="Path to tab-separated summary stats")
parser$add_argument("--bim", type="character", required=TRUE, help="Path to genotype .bim file for allele alignment")
parser$add_argument("--output", type="character", required=TRUE, help="Path for the output file")
parser$add_argument("--n_total", type="integer", help="Total sample size (if constant)")

args <- parser$parse_args()

# 2. LOAD DATA
message("Reading input file: ", args$input)
# fread2 handles tab-separated files automatically
df <- bigreadr::fread2(args$input)

message("Reading BIM file: ", args$bim)
bim <- bigreadr::fread2(args$bim, select = c(1, 2, 4, 5, 6))
colnames(bim) <- c("bim.chr", "rsid", "bim.pos", "bim.a1", "bim.a0")

# 3. TRANSFORMATION LOGIC
message("Aligning alleles by RSID and calculating n_eff if not provided...")

processed_df <- df %>%
  # Match your specific column names to internal variables
  rename(rsid = rsids) %>%
  # Join with BIM to ensure genomic coordinates match your genotypes
  inner_join(bim, by = c("rsid")) %>%
  mutate(
    n_eff = if (!is.null(args$n_total)) {
      args$n_total
    } else {
      1 / (2 * af * (1 - af) * sebeta^2)
    },

    # Standardize Reference Allele (a0)
    a1 = alt,
    a0 = ifelse(a1 == bim.a1, bim.a0, bim.a1),

    # Ensure P-value is present (using your pval column)
    CHR = as.integer(bim.chr),
    BP = as.integer(bim.pos),
    P = pval
  ) %>%
  # Format specifically for PRS tools
  select(
    SNP = rsid,
    CHR = CHR,
    BP = BP,
    A1 = a1,
    A2 = a0,
    BETA = beta,
    SE = sebeta,
    P = P,
    n_eff = n_eff
  ) %>%
  # Remove any rows where n_eff calculation resulted in Inf or NA due to AF=0 or 1
  filter(is.finite(n_eff))

# 4. SAVE OUTPUT
message("Writing cleaned stats to: ", args$output)
write.table(processed_df, 
            file = args$output, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

message("Done! Processed ", nrow(processed_df), " SNPs.")