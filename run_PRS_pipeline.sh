#!/bin/bash 
#SBATCH --time=13:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH -o prs_pipeline.out
#SBATCH --job-name prs_pipeline


# ---- Defaults (documented only; must be overridden or validated) ----
# Needs the conda environment gdcPipeline
path_code="/projects/standard/gdc/public/prs_methods/scripts/PRScsx"
path_ref_dir="/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref"
anc1="AFR"
anc2="EUR"
target_sumstats_file="${path_data_root}/gwas/target_sumstats.txt"
training_sumstats_file="${path_data_root}/gwas/training_sumstats.txt"
reference_SNPS_bim="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample" # This is the full path except for extension for bim file
study_sample_plink_anc2="${path_data_root}/anc2_plink_files/${anc2}_simulation_study_sample"
prs_pipeline="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline" # Path to the prs_pipeline code

####### Variables ####### # VIPRS
# Needs the conda environment viprs_env
bfile_gwas_input=${path_data}/anc1_plink_files/archived/AFR_simulation_gwas
covariate_file_gwas=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/prs_pipeline/viprs/gwas/temp/viprs_summary_stats_covar_sex_no_header.txt
covariate_file_study_sample=${path_data}/prs_pipeline/viprs/study_sample_covar.txt # No header FID IID Sex

##### Common variables #####
path_plink2=/projects/standard/gdc/public/plink2
path_data=/projects/standard/gdc/public/prs_methods/data/simulated_1000G # AKA path_data_root (prscsx)
out_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G # AKA output_dir (prscsx)
study_sample_plink="${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample" # AKA bfile_study_sample (prscsx)

# to make cross-compatable
path_data_root="${path_data}"
output_dir="${out_path}"
bfile_study_sample="${study_sample_plink}"
