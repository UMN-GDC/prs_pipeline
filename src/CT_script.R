#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4){
  stop("Usage: run_prs.R <phenotype> <pcs> <profile_prefix> <final_output_dir>")
}

phenotype_file <- args[1]
pcs_file <- args[2]
profile_prefix <- args[3]
final_output_dir <- args[4]

p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)

# Read phenotype
phenotype <- read.table(phenotype_file, header=TRUE)

# Read PCs
pcs <- read.table(pcs_file, header=FALSE)
colnames(pcs) <- c("FID","IID", paste0("PC", 1:(ncol(pcs)-2)))

# Merge phenotype + PCs
pheno <- merge(phenotype, pcs, by=c("FID","IID"))

# Null model (no PRS)
null.model <- lm(phenotype~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
null.r2 <- summary(null.model)$r.squared

prs.result <- NULL

for(i in p.threshold){

  profile_file <- paste0(profile_prefix, ".", i, ".profile")

  prs <- read.table(profile_file, header=TRUE)

  pheno.prs <- merge(pheno, prs[,c("FID","IID","SCORE")], by=c("FID","IID"))

  model <- lm(phenotype~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])

  model.r2 <- summary(model)$r.squared
  prs.r2 <- model.r2 - null.r2

  prs.coef <- summary(model)$coeff["SCORE",]

  prs.beta <- as.numeric(prs.coef[1])
  prs.se <- as.numeric(prs.coef[2])
  prs.p <- as.numeric(prs.coef[4])

  prs.result <- rbind(
    prs.result,
    data.frame(
      Threshold=i,
      R2=prs.r2,
      P=prs.p,
      BETA=prs.beta,
      SE=prs.se
    )
  )
}

dir.create(final_output_dir, showWarnings=FALSE, recursive=TRUE)

CT.output <- file.path(final_output_dir, "CT_prs_results.txt")

write.table(prs.result,
            file=CT.output,
            sep=" ",
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE)

# Print best result to stdout
print(prs.result[which.max(prs.result$R2),])
