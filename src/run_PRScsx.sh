#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=5:00:00
#SBATCH -p msismall
#SBATCH -o PRScsx_sims.out
#SBATCH --job-name PRScsx_sims
# #SBATCH --mail-type=ALL  
# #SBATCH --mail-user=baron063@umn.edu 


#### Load environment ####
source /projects/standard/gdc/public/envs/load_miniconda3.sh

#### Initial script creation to be modified to be command line callable.
#### Needs to have the flexibility to help combine the anc1 and anc2 plink files into a combined dataset if user desires


#### Variables ####
path_code=/projects/standard/gdc/public/prs_methods/scripts/PRScsx # Path to where PRScsx is cloned
path_data_root=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_ref_dir=/projects/standard/gdc/public/prs_methods/ref/ref_PRScsx/1kg_ref
path_plink2=/projects/standard/gdc/public/plink2

anc1=AFR # This turns into your target ancestry... 
anc2=EUR # This turns into your training ancestry...
nsims=1 # This is the number of simulations you'd like generated. 

target_sumstats_file=${path_data_root}/gwas/target_sumstats.txt    # rsid    chr     a1      a0      beta    beta_se n_eff
training_sumstats_file=${path_data_root}/gwas/training_sumstats.txt # SNP A1 A2 BETA SE (tab separated with a header line)

output_dir=${path_data_root}/prs_pipeline/PRScsx

## Benefits from a combined train_validate data split
validation_bim=${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample_train_validate # full path and prefix for the target dataset 
test_bim=${path_data_root}/anc1_plink_files/${anc1}_simulation_study_sample_test



#### Derived variables ####
## Renaming columns for PRScsx useage ##
target_sst_file_using=${path_data_root}/gwas/target_sumstats_PRScsx.txt # SNP A1 A2 BETA SE (tab separated with a header line)

awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","A1","A2","BETA","SE"} NR>1 {print $1,$3,$4,$5,$6}' "${target_sumstats_file}" > "${target_sst_file_using}"

## Extract from tlprs sumstats _ sim num the N value for below?
GWAS_sample_size_target=$(awk 'NR==2 {print $NF}' "$target_sumstats_file")

## For training sumstats ##
training_sst_file_using=${path_data_root}/gwas/training_sumstats_PRScsx.txt # SNP A1 A2 BETA SE (tab separated with a header line)

awk 'BEGIN {OFS="\t"} NR==1 {print "SNP","A1","A2","BETA","SE"} NR>1 {print $1,$3,$4,$5,$6}' "${training_sumstats_file}" > "${training_sst_file_using}"

## Extract from tlprs sumstats _ sim num the N value for below?
GWAS_sample_size_training=$(awk 'NR==2 {print $NF}' "$training_sumstats_file")

#### Script start ####
mkdir -p ${output_dir}

echo "Python call:"
echo "python ${path_code}/PRScsx.py --ref_dir="${path_ref_dir}" \
  --bim_prefix=${validation_bim} \
  --sst_file=${training_sst_file_using},${target_sst_file_using} \
  --n_gwas=${GWAS_sample_size_training},${GWAS_sample_size_target} \
  --pop=${anc2},${anc1} \
  --out_dir=${output_dir} \
  --out_name=PRScsx_joint_${anc1} \
  --seed=42"

python ${path_code}/PRScsx.py --ref_dir="${path_ref_dir}" \
  --bim_prefix=${validation_bim} \
  --sst_file=${training_sst_file_using},${target_sst_file_using} \
  --n_gwas=${GWAS_sample_size_training},${GWAS_sample_size_target} \
  --pop=${anc2},${anc1} \
  --out_dir=${output_dir} \
  --out_name=PRScsx_joint_${anc1} \
  --seed=42

## Combine per-chromosome PRScsx output into a single file ##
pushd ${output_dir}
out_prefix="PRScsx_joint_${anc1}_${anc1}_pst_eff_a1_b0.5_phiauto"
combined_file="PRScsx_joint_${anc1}_${anc1}_weights.txt"

# Create header
echo -e "SNP\tA1\tBETA" > $combined_file

# Append data from each chromosome (skip headers) #Should have all 6 columns
for chr in {1..22}; do
  awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6}' ${out_prefix}_chr${chr}.txt >> $combined_file
done


## Generate scores ##
${path_plink2} --bfile ${test_bim} --score ${combined_file} 2 4 6 header --out PRScsx_joint_${anc1}_score


#### Calculate R squared value ####
Rscript /projects/standard/gdc/shared/code/TLPRS/TLPRS_tools/src/PRS_sscore_to_R2.R PRScsx_joint_${anc1}_score.sscore
cat adj_PRS_sscore_Rsqr.txt
popd

