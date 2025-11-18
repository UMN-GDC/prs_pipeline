#!/usr/bin/env Rscript

# ==========================================================
# Title: Split Target Population Data for TLPRS (Extended Version)
# Author: Kody DeGolier
# Description:
#   Splits a target population PLINK dataset (.fam/.bed/.bim)
#   and corresponding simulated phenotype data into GWAS,
#   training, validation, and testing subsets.
#
# Inputs:
#   1. location          - Directory containing data files
#   2. filename_1        - Target population .fam file
#   3. output_name       - (Optional) Base output name (no extension)
#   4. gwas_num          - (Optional) Proportion for GWAS (0-1)
#   5. train_num         - (Optional) Proportion for training (0-1)
#   6. validate_num      - (Optional) Proportion for validation (0-1)
#
# Outputs:
#   - *_gwasdata.csv, *_traindata.csv, *_validatedata.csv, *_testdata.csv
#   - *_simulated.txt, *_simulated_gwas.csv, *_simulated_train.csv, *_simulated_validate.csv, *_simulated_test.csv
#   - PLINK datasets for each split (_gwas_split, _training_split, _validate_split, _testing_split)
#
# Example:
#   Rscript split_target_data.R /path/to/data target.fam target 0.4 0.3 0.1
# ==========================================================

# ----------------------------------------------------------
# Load Packages
# ----------------------------------------------------------
.libPaths("/home/gdc/public/Ref/R")
setwd("/home/gdc/public/Ref/R")
suppressMessages(library(tidyverse))

# ----------------------------------------------------------
# Parse Command-Line Arguments
# ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript split_target_data.R <location> <filename_1> [output_name] [gwas_num] [train_num] [validate_num]", call. = FALSE)
}

location  <- args[1]
filename_1 <- args[2]
base_filename_1 <- str_sub(filename_1, 1, -5) # remove .fam extension

# Default parameters if not provided
if (length(args) == 2) {
  output_name <- "default_name"
  gwas_percent <- 0.4
  train_percent <- 0.3
  validate_percent <- 0.1
} else {
  output_name <- args[3]
  gwas_percent <- as.numeric(args[4])
  train_percent <- as.numeric(args[5])
  validate_percent <- as.numeric(args[6])
}

test_percent <- 1 - gwas_percent - train_percent - validate_percent
if (test_percent <= 0) stop("Proportions must sum to less than 1.")

# ----------------------------------------------------------
# Helper Function for Logging
# ----------------------------------------------------------
log_msg <- function(...) cat("[INFO]", ..., "\n")

# ----------------------------------------------------------
# Initialization
# ----------------------------------------------------------
setwd(location)
log_msg("Working directory:", location)
log_msg("Input .fam file:", filename_1)
log_msg("Output base name:", output_name)
log_msg(sprintf("Split proportions  GWAS: %.2f | Train: %.2f | Validate: %.2f | Test: %.2f",
                gwas_percent, train_percent, validate_percent, test_percent))

# ----------------------------------------------------------
# Step 1: Read Input .fam File
# ----------------------------------------------------------
if (!file.exists(filename_1)) stop(paste("Input .fam file not found:", filename_1))

file1_dat <- read.table(filename_1, header = FALSE)
colnames(file1_dat) <- c("FID", "IID", "Fatherinstdy", "Motherinstdy", "Sex", "Y")

# ----------------------------------------------------------
# Step 2: Split by Individual IDs
# ----------------------------------------------------------
set.seed(10)
all_ids <- file1_dat$IID

# Sample GWAS IDs first
gwas_ids <- sample(all_ids, size = gwas_percent * length(all_ids), replace = FALSE)
remaining_1 <- file1_dat %>% filter(!IID %in% gwas_ids)

# Sample training IDs
train_ids <- sample(remaining_1$IID, size = train_percent * length(all_ids), replace = FALSE)
remaining_2 <- remaining_1 %>% filter(!IID %in% train_ids)

# Sample validation IDs
validate_ids <- sample(remaining_2$IID, size = validate_percent * length(all_ids), replace = FALSE)
test_ids <- setdiff(remaining_2$IID, validate_ids)

# Create subset dataframes
gwas_data <- file1_dat %>% filter(IID %in% gwas_ids) %>% select(FID, IID)
train_data <- file1_dat %>% filter(IID %in% train_ids) %>% select(FID, IID)
validate_data <- file1_dat %>% filter(IID %in% validate_ids) %>% select(FID, IID)
test_data <- file1_dat %>% filter(IID %in% test_ids) %>% select(FID, IID)

# ----------------------------------------------------------
# Step 3: Write Split ID Lists
# ----------------------------------------------------------
write.table(gwas_data, paste0(output_name, "_gwasdata.csv"), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(train_data, paste0(output_name, "_traindata.csv"), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(validate_data, paste0(output_name, "_validatedata.csv"), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(test_data, paste0(output_name, "_testdata.csv"), sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

log_msg("Wrote PLINK ID splits (.csv files)")

# ----------------------------------------------------------
# Step 4: Process Phenotype Data
# ----------------------------------------------------------
phenotype_path <- paste0(base_filename_1, ".fam")
if (!file.exists(phenotype_path)) stop(paste("Phenotype .fam file not found:", phenotype_path))

phenotype <- read.table(phenotype_path, header = FALSE)
colnames(phenotype) <- c("FID", "IID", "V3", "V4", "Sex", "Y")
phenotype <- phenotype %>% select(-c(V3, V4))

# Write full simulated phenotype
write.table(phenotype, paste0(output_name, "_simulated.txt"), sep = " ", quote = FALSE, row.names = FALSE)

# Write phenotype splits
write.table(phenotype %>% filter(IID %in% gwas_ids),
            paste0(output_name, "_simulated_gwas.csv"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(phenotype %>% filter(IID %in% train_ids),
            paste0(output_name, "_simulated_train.csv"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(phenotype %>% filter(IID %in% validate_ids),
            paste0(output_name, "_simulated_validate.csv"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(phenotype %>% filter(IID %in% test_ids),
            paste0(output_name, "_simulated_test.csv"), sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)

log_msg("Wrote phenotype splits")

# ----------------------------------------------------------
# Step 5: Run PLINK Commands
# ----------------------------------------------------------
plink_base_cmd <- "plink --bfile "
plink_keep <- " --keep "
plink_suffix <- " --make-bed --out "

gwas_out <- paste0(output_name, "_gwas_split")
train_out <- paste0(output_name, "_training_split")
val_out   <- paste0(output_name, "_validate_split")
test_out  <- paste0(output_name, "_testing_split")

plink_call1 <- paste0(plink_base_cmd, base_filename_1, plink_keep, output_name, "_gwasdata.csv", plink_suffix, gwas_out)
plink_call2 <- paste0(plink_base_cmd, base_filename_1, plink_keep, output_name, "_traindata.csv", plink_suffix, train_out)
plink_call3 <- paste0(plink_base_cmd, base_filename_1, plink_keep, output_name, "_validatedata.csv", plink_suffix, val_out)
plink_call4 <- paste0(plink_base_cmd, base_filename_1, plink_keep, output_name, "_testdata.csv", plink_suffix, test_out)

log_msg("Running PLINK splits...")
system(plink_call1)
system(plink_call2)
system(plink_call3)
system(plink_call4)

# ----------------------------------------------------------
# Completion
# ----------------------------------------------------------
log_msg("Data successfully split into GWAS, training, validation, and testing sets.")
log_msg("All PLINK-ready files generated.")
