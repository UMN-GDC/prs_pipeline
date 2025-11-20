#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=50g
#SBATCH --tmp=50g
#SBATCH -p agsmall
#SBATCH -o generate_summary_stats.out
#SBATCH --job-name prs_generate_summary_stats

# set -euo pipefail
# IFS=$'\n\t'

# ========== CONFIGURATION ==========
base_location="/home/gdc/public/prs_methods/data/test/sim_1"
anc1_gwas_input="AFR_simulation_gwas"
anc2_gwas_input="EUR_simulation_gwas"

base_location="/home/gdc/public/prs_methods/outputs/simulation_AFR_target_EUR_training_simPheno_test/sim_1"

path_to_repo="/home/gdc/public/prs_methods/scripts/prs_pipeline"

log_file="${base_location}/gen_sum_stats.log"

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


# ========== STEP 1: Generate .qassoc and summary stats ==========
log "STEP 1: Generating association statistics for provided gwas_split data"

plink --bfile "${base_location}/${anc1_gwas_input}" --assoc --allow-no-sex \
  --out "${base_location}/ancestry_target" \
  &>> "${log_file}"

plink --bfile "${base_location}/${anc2_gwas_input}" --assoc --allow-no-sex \
  --out "${base_location}/ancestry_training" \
  &>> "${log_file}"


# ========== STEP 2: PCA ==========
log "STEP 2: Performing PCA for ${anc1_gwas_input} and ${anc2_gwas_input}"

plink --bfile "${base_location}/${anc1_gwas_input}" --pca --allow-no-sex \
  --out "${base_location}/${anc1_gwas_input}_pca" &>> "${log_file}"

plink --bfile "${base_location}/${anc2_gwas_input}" --pca --allow-no-sex \
  --out "${base_location}/${anc2_gwas_input}_pca" &>> "${log_file}"

# ========== STEP 3: GWAS summary stats corrected ==========
log "STEP 3: Running GWAS summary stats with PCA covariates"

plink --bfile "${base_location}/${anc1_gwas_input}" \
  --covar "${base_location}/${anc1_gwas_input}_pca.eigenvec" \
  --linear --allow-no-sex \
  --out "${base_location}/target_sumstats_corrected" &>> "${log_file}"

plink --bfile "${base_location}/${anc2_gwas_input}" \
  --covar "${base_location}/${anc2_gwas_input}_pca.eigenvec" \
  --linear --allow-no-sex \
  --out "${base_location}/training_sumstats_corrected" &>> "${log_file}"

# ========== STEP 4: Filter and format summary stats ==========
log "STEP 4: Filtering and formatting summary stats"

for type in training target; do
  input="${base_location}/${type}_sumstats_corrected.assoc.linear"
  filtered="${base_location}/${type}_sumstats_filtered.txt"
  clean="${base_location}/${type}_sumstats_filtered_clean.txt"

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
  "target_sumstats_filtered_clean.txt" \
  "${base_location}/${anc1_gwas_input}.bim" \
  "target_sumstats.txt" \
  "training_sumstats_filtered_clean.txt" \
  "${base_location}/${anc2_gwas_input}.bim" \
  "training_sumstats.txt" \
  > "${base_location}/creating_target_sumstats_files_for_TLPRS.log" 2>&1


Rscript "${path_to_repo}/src/creating_target_sumstats_files_multiple.R" \
  "${base_location}" \
  "target_sumstats_filtered_clean.txt" \
  "${base_location}/${anc1_gwas_input}.bim" \
  "target_sumstats_prosper.txt" \
  "training_sumstats_filtered_clean.txt" \
  "${base_location}/${anc2_gwas_input}.bim" \
  "training_sumstats_prosper.txt" \
  > "${base_location}/creating_target_sumstats_files.log" 2>&1


log "Finished generating batch of summary statistic files"
