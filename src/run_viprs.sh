#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=48:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=baron063@umn.edu 
#SBATCH -o VIprs_sims.out
#SBATCH --job-name VIprs_sims

source /projects/standard/gdc/public/envs/load_miniconda3.sh 
conda deactivate
conda activate viprs_env

module load plink


####### CONFIG #######

path_data=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
targ_ancestry=AFR
out_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_prs_pipeline=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline

# Inputs
piece1=${path_data}/gwas/target_sumstats_CTSLEB.txt      # CHR SNP BP A1 beta beta_se p rs_id
piece2=${path_data}/gwas/target_sumstats.txt             # rsid chr a1 a0 beta beta_se n_eff
piece3=${path_data}/gwas/target_sumstats_corrected.assoc.linear # CHR SNP BP A1 TEST NMISS BETA STAT P

# Output
out_file=${path_data}/prs_pipeline/viprs/${targ_ancestry}_sumstats_viprs.txt

source ${path_prs_pipeline}/src/viprs_helper.sh

make_viprs_sumstats ${piece1} ${piece2} ${piece3} ${out_path}/prs_pipeline/viprs/gwas

viprs_fit -l "/projects/standard/gdc/public/prs_methods/ref/ref_viprs/EUR/chr_*" \
  -s ${piece2} \
  --sumstats-format custom \
  --custom-sumstats-mapper \
    rsid=rsid,chr=chr,eff_allele=a1,non_eff_allele=a0,beta=beta,se=best_se,n_eff=n_eff \
  --output-dir "${out_path}/prs_pipeline/viprs" \
  --model VIPRS \
  --exclude-lrld \
  --float-precision float64

# for chr in {1..22}; do
#     sed -i '1s/^#//' /projects/standard/gdc/public/prs_methods/outputs/simulation_AFR_target_EUR_training_simPheno/viprs_by_chr/chr${chr}.txt
# done


python stable_archive/viprs_evaluate_pseudo_r2.py \
  --sumstats_files "${provided_sumstats}" \
  --sumstats_format plink \
  --inferred_betas "/projects/standard/gdc/public/prs_methods/outputs/viprs/VIPRS_EM.fit.gz"
