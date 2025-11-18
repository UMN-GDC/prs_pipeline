#!/usr/bin/env bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=8
#SBATCH --mem=16g
#SBATCH --tmp=32g
#SBATCH -p agsmall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=baron063@umn.edu 
#SBATCH -o run_simulation_generate_inputs.out
#SBATCH --job-name prs_generate_summary_stats

# set -euo pipefail
# IFS=$'\n\t'

# ========== CONFIGURATION ==========
counter=1
anc1="AFR"
anc2="EUR"
output_name_target="AFR_simulation_output"
match_prefix_train_pop="EUR_sim_pheno.QC8_"
path_to_repo="/home/gdc/shared/code/TLPRS/TLPRS_tools"
base_location="/home/gdc/public/prs_methods/outputs/simulation_AFR_target_EUR_training_simPheno"
output_prefix_anc1="/${base_location}/${anc1}_simulation_output_gwas_split"
output_prefix_anc2="/${base_location}/${anc2}_sim_pheno.QC8_gwas"

log_file="${base_location}/pipeline_run_${counter}.log"

if [ -f ${log_file} ]; then
  rm ${log_file}
fi

mkdir -p "${base_location}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${log_file}"
}

# ========== Modules ==========
module load R/4.4.0-openblas-rocky8
module load plink
module load plink/2.00-alpha-091019

success=0

# ========== STEP 1: Generate .qassoc and summary stats ==========
log "STEP 1: Generating association statistics for run ${counter}"

plink --bfile "${base_location}/${output_name_target}_gwas_split" --assoc --allow-no-sex \
  --out "${base_location}/target_${counter}" \
  &>> "${log_file}"

plink --bfile "${base_location}/${match_prefix_train_pop}gwas" --assoc --allow-no-sex \
  --out "${base_location}/training_${counter}" \
  &>> "${log_file}"

# Rename to standard names
mv "${base_location}/target_${counter}.qassoc" "${base_location}/ancestry_target_${counter}.qassoc"
mv "${base_location}/training_${counter}.qassoc" "${base_location}/ancestry_training_${counter}.qassoc"

# ========== STEP 2: PCA ==========
log "STEP 2: Performing PCA for ${anc1} and ${anc2}"

plink --bfile "${output_prefix_anc1}" --pca --allow-no-sex \
  --out "${base_location}/${anc1}_QC8_pca_${counter}" &>> "${log_file}"

plink --bfile "${output_prefix_anc2}" --pca --allow-no-sex \
  --out "${base_location}/${anc2}_QC8_pca_${counter}" &>> "${log_file}"

# ========== STEP 3: GWAS summary stats corrected ==========
log "STEP 3: Running GWAS summary stats with PCA covariates"

plink --bfile "${output_prefix_anc1}" \
  --covar "${base_location}/${anc1}_QC8_pca_${counter}.eigenvec" \
  --linear --allow-no-sex \
  --out "${base_location}/target_sumstats_corrected_${counter}" &>> "${log_file}"

plink --bfile "${output_prefix_anc2}" \
  --covar "${base_location}/${anc2}_QC8_pca_${counter}.eigenvec" \
  --linear --allow-no-sex \
  --out "${base_location}/training_sumstats_corrected_${counter}" &>> "${log_file}"

# ========== STEP 4: Filter and format summary stats ==========
log "STEP 4: Filtering and formatting summary stats"

for type in training target; do
  input="${base_location}/${type}_sumstats_corrected_${counter}.assoc.linear"
  filtered="${base_location}/${type}_sumstats_filtered_${counter}.txt"
  clean="${base_location}/${type}_sumstats_filtered_clean_${counter}.txt"

  awk 'BEGIN {OFS="\t"}
    NR==1 { print "SNP","CHR","A1","BETA","SE","P","n_eff"; next }
    $5=="ADD" {
      se = ($8==0 || $8=="NA") ? "NA" : $7/$8
      print $2,$1,$4,$7,se,$9,$6
    }' "${input}" > "${filtered}"

  awk '$5 != "NA"' "${filtered}" > "${clean}"
done

# ========== STEP 5: Create target summary stats files ==========
log "STEP 5: Creating target summary stats files via R"

Rscript "${path_to_repo}/src/creating_target_sumstats_files.R" \
  "${base_location}" \
  "target_sumstats_filtered_clean_${counter}.txt" \
  "${base_location}/${output_name_target}_gwas_split.bim" \
  "target_sumstats_${counter}.txt" \
  "training_sumstats_filtered_clean_${counter}.txt" \
  "${base_location}/${match_prefix_train_pop}gwas.bim" \
  "training_sumstats_${counter}.txt" \
  > "${base_location}/creating_target_sumstats_files_${counter}_for_TLPRS.log" 2>&1

Rscript "${path_to_repo}/src/creating_target_sumstats_files_multiple.R" \
  "${base_location}" \
  "target_sumstats_filtered_clean_${counter}.txt" \
  "${base_location}/${output_name_target}_gwas_split.bim" \
  "target_sumstats_prosper_${counter}.txt" \
  "training_sumstats_filtered_clean_${counter}.txt" \
  "${base_location}/${match_prefix_train_pop}gwas.bim" \
  "training_sumstats_prosper_${counter}.txt" \
  > "${base_location}/creating_target_sumstats_files_${counter}.log" 2>&1


success=$((success + 1))
log "Run ${counter} completed successfully. Success count: ${success}"
