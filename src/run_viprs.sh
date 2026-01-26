#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=12:00:00
#SBATCH -p msismall
#SBATCH -o VIprs_sims.out
#SBATCH --job-name VIprs_sims

source /projects/standard/gdc/public/envs/load_miniconda3.sh 
conda deactivate
conda activate viprs_env

module load plink


###### FUNCTION ######
generate_viprs_sumstats() {
local path_plink2="$1"
local bfile_input="$2"
local out_name="$3"
local covariate_file="$4"

for chr in {1..22}; do
  ${path_plink2} \
    --bfile "${bfile_input}" \
    --covar "${covariate_file}" \
    --chr $chr \
    --glm hide-covar cols=+a1freq,+nobs \
    --out "${out_name}"_${chr}

  awk '
  BEGIN {
    OFS="\t";
    print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
  }
  NR>1 && $10=="ADD" {
    print $1,$2,$3,$4,$5,$7,$9,$11,$12,$13,$14,$15
  }
  ' "${out_name}_${chr}.PHENO1.glm.linear" > "${out_name}_${chr}.PHENO1_corrected.glm.linear"
done
}

load_config() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        echo "Loading configuration from: $file_path"
        # Source the file. Variables set in the config will override defaults.
        # Note: 'source' is used for shell-readable KEY="VALUE" files.
        source "$file_path"
        return 0
    else
        echo "ERROR: Configuration file not found at $file_path" >&2
        return 1
    fi
}

####### Variables #######
path_data=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
targ_ancestry=AFR
out_path=/projects/standard/gdc/public/prs_methods/data/simulated_1000G
path_plink2=/projects/standard/gdc/public/plink2
bfile_gwas_input=${path_data}/anc1_plink_files/archived/AFR_simulation_gwas
bfile_study_sample=${path_data}/anc1_plink_files/AFR_simulation_study_sample
covariate_file_gwas=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/prs_pipeline/viprs/gwas/temp/viprs_summary_stats_covar_sex_no_header.txt
covariate_file_study_sample=${path_data}/prs_pipeline/viprs/study_sample_covar.txt

# Output
out_file=${out_path}/prs_pipeline/viprs/viprs_summary_stats
out_dir=$(dirname "$out_file")


## If the user doesn't provide a covariate file
# awk 'BEGIN{OFS="\t"; print "FID","IID","SEX"} {print $1,$2,$5}' ${bfile_gwas_input}.fam > ${out_file}_covar_sex.txt

#### # ---- Parse args ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --c) config_file="$2"; shift 2; break ;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

if [[ -n "$config_file" ]]; then
    # Case 1: --config was specified. Load it, and it overrides ALL defaults.
    load_config "$config_file" || exit 1
else
    # Case 2: No --config was specified.
    # The variables set by the individual flags in step 3 are used.
    echo "No --c file specified. Using command-line arguments and defaults."
fi


generate_viprs_sumstats ${path_plink2} ${bfile_gwas_input} ${out_file} ${covariate_file_gwas}

mkdir -p ${out_dir}/gwas/logs
mv ${out_dir}/*.log ${out_dir}/gwas/logs
mv ${out_dir}/*PHENO* ${out_dir}/gwas
mkdir -p ${out_dir}/gwas/temp
mv ${out_dir}/gwas/*PHENO1.* ${out_dir}/gwas/temp

provided_sumstats="${out_dir}/gwas/viprs_summary_stats_*"

viprs_fit -l "/projects/standard/gdc/public/prs_methods/ref/ref_viprs/AFR/chr_*" \
  -s "${provided_sumstats}" \
  --sumstats-format plink \
  --output-dir "${out_dir}" \
  --model VIPRS \
  --exclude-lrld \
  --float-precision float64


echo "Preview of outputed posterior distribution of the variant effect sizes produced from viprs_fit call:"
zcat ${out_dir}/VIPRS_EM.fit.gz | head

viprs_score -f "${out_dir}/VIPRS_EM.fit.gz" \
             --bfile "${bfile_study_sample}" \
             --output-file "${out_dir}/VIPRS_PGS"

echo "Preview of outputed PRS produced from viprs_score call:"
head ${out_dir}/VIPRS_PGS.prs

#### Making the pheno file assuming it is not provided ####
awk 'BEGIN{OFS="\t"} {print $1,$2,$6}' ${bfile_study_sample}.fam > ${out_dir}/study_sample_pheno.txt

# #### Making the testing covariate file if not provided ####
#awk 'BEGIN{OFS="\t"} {print $1,$2,$5}' ${bfile_study_sample}.fam > ${out_dir}/study_sample_covar.txt

viprs_evaluate --prs-file "${out_dir}/VIPRS_PGS.prs" \
                --phenotype-file "${out_dir}/study_sample_pheno.txt" \
                --phenotype-col 2 \
                --covariates-file "${covariate_file_study_sample}" \
                --output-file ${out_dir}/viprs_evaluate_results

echo "Checking results from viprs_evaluate on the test data"
head ${out_dir}/viprs_evaluate_results.eval
