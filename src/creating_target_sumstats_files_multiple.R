#!/usr/bin/env Rscript
#Obtaining necessary package
.libPaths("/home/gdc/public/Ref/R")
suppressMessages(library(readr))
suppressMessages(library(tidyverse))


# --- helpers ---------------------------------------------------------------
read_assoc <- function(path) {
  read_table(path, col_types = cols()) %>%
    transmute(
      rsid = SNP,
      chr = CHR,
      effect_allele_assoc = A1,
      beta = BETA,
      beta_se = SE,
      p = P,
      n_eff = n_eff
    )
}

read_bim <- function(path) {
  read_table(path, col_names = FALSE, col_types = cols()) %>%
    # PLINK .bim: chr, rsid, cm, bp, A1, A2
    set_names(c("chrom", "rsid", "cm", "bp", "a1_bim", "a0_bim")) %>%
    transmute(
      rsid = rsid,
      a1_bim = a1_bim,
      a0_bim = a0_bim,
      bp = bp
    )
}

# Merge assoc + bim, flipping beta if effect allele disagrees
merge_with_bim <- function(assoc_df, bim_df, keep_a0 = TRUE) {
  assoc_df %>%
    inner_join(bim_df, by = "rsid") %>%
    mutate(
      flip = effect_allele_assoc != a1_bim,
      beta = ifelse(flip, -beta, beta),
      a1 = effect_allele_assoc,  # after resolution, treat assoc A1 as effect allele
      a0 = ifelse(flip, a1_bim, a0_bim)  # keep the non-effect allele consistently
    ) %>%
    {
      if (keep_a0) {
        select(., rsid, chr, a1, a0, beta, beta_se, n_eff)
      } else {
        select(., rsid, chr, a1, beta, beta_se, n_eff)
      }
    } %>%
    na.omit()
}

write_ctsleb <- function(assoc_df, bim_df, out_path) {
  # Need CHR, SNP, BP, A1, BETA, SE, P, rs_id
  assoc_df %>%
    inner_join(bim_df, by = "rsid") %>%
    mutate(
      CHR = chr,
      SNP = rsid,
      BP = bp,
      A1 = effect_allele_assoc,
      rs_id = rsid,
      SNP_2 = paste(rsid, BP, A1, sep = ":")
    ) %>%
    select(CHR, SNP, BP, A1, beta, beta_se, p, rs_id) %>%
    na.omit() %>%
    write.table(file = out_path, sep = " ", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
}

# --- argument parsing -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(4, 7))) {
  stop("This script expects either 4 or 7 arguments.", call. = FALSE)
}

base_location          <- args[1] # Path to data
simulated_assoc_name_1 <- args[2] # Target summarystats base
simulated_bim_name_1   <- args[3] # Target .bim file
sum_stat_1_name        <- args[4] # Output file for target summarystats

has_training <- length(args) == 7
if (has_training) {
  simulated_assoc_name_2 <- args[5] # Training summarystats base
  simulated_bim_name_2   <- args[6] # Training .bim file
  sum_stat_2_name        <- args[7] # Output file for training summarystats
}

setwd(base_location)

# --- target / "prosper" sumstats ------------------------------------------
assoc1 <- read_assoc(simulated_assoc_name_1)
bim1  <- read_bim(simulated_bim_name_1)

prosper_sumstats <- merge_with_bim(assoc1, bim1, keep_a0 = TRUE)
write.table(prosper_sumstats, file = sum_stat_1_name,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# --- CTSLEB target --------------------------------------------------------
base_ss_name <- tools::file_path_sans_ext(sum_stat_1_name)
ctsleb_targ_name <- paste0(base_ss_name, "_CTSLEB.txt")
write_ctsleb(assoc1, bim1, ctsleb_targ_name)

# --- optional training ----------------------------------------------------
if (has_training) {
  assoc2 <- read_assoc(simulated_assoc_name_2)
  bim2  <- read_bim(simulated_bim_name_2)

  training_sumstats <- merge_with_bim(assoc2, bim2, keep_a0 = TRUE)
  write.table(training_sumstats, file = sum_stat_2_name,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  # CTSLEB training
  base_ss_name_2 <- tools::file_path_sans_ext(sum_stat_2_name)
  ctsleb_train_name <- paste0(base_ss_name_2, "_CTSLEB.txt")
  write_ctsleb(assoc2, bim2, ctsleb_train_name)
}
