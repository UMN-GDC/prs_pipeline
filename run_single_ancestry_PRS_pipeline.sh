#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH -p agsmall
#SBATCH -o prs_unified_%j.out
#SBATCH --job-name prs_pipeline


# --- Defaults ---
RUN_CT=false
RUN_LDPRED2=false
RUN_LASSOSUM2=false

# Paths (Consolidated)
summary_stats_file="/projects/standard/gdc/public/prs_methods/data/Phenotypes/updated_summary_stats/updated_phenocode-aseg_lh_volume_Hippocampus.tsv"
bim_file_path="/projects/standard/gdc/public/prs_methods/data/simulated_1000G/anc1_plink_files/archived/AFR_simulation_gwas.bim"
study_sample="/projects/standard/gdc/public/prs_methods/data/simulated_1000G/anc1_plink_files/AFR_simulation_study_sample"
output_path="/projects/standard/gdc/public/prs_methods/data/simulated_1000G"
path_repo="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
n_total_gwas=33224
gwas_pca_eigenvec_file="/projects/standard/gdc/shared/abcdTest/04-globalAncestry/merged_dataset_pca.eigenvec"

# --- Environment ---
module load R/4.4.0-openblas-rocky8
export R_LIBS_USER="/projects/standard/gdc/public/Ref/R"

# --- Pre-load config ---
for arg in "$@"; do
    if [[ "$arg" == "-C" ]]; then
        # This gets the value immediately following -C
        conf_file=$(echo "$@" | grep -oP '(?<=-C )[^ ]+')
        if [[ -f "$conf_file" ]]; then
            echo "Loading config: $conf_file"
            source "$conf_file"
        else
            echo "Error: Config file $conf_file not found."
            exit 1
        fi
    fi
done

# --- Usage ---
usage() {
    echo "Usage: $0 [-C config_file] [-c] [-l] [-s]"
    echo "  -C    Path to config file"
    echo "  -c    Run Clumping + Thresholding (C+T)"
    echo "  -l    Run LDpred2"
    echo "  -s    Run lassosum2"
    exit 1
}

# --- Arg Parser ---
while getopts "C:cls" opt; do
    case "$opt" in
        C) ;; # Handled in pre-load step, but kept here so getopts doesn't complain
        c) RUN_CT=true ;;
        l) RUN_LDPRED2=true ;;
        s) RUN_LASSOSUM2=true ;;
        *) usage ;;
    esac
done

# If no flags provided, show usage
if [[ "$RUN_CT" == false && "$RUN_LDPRED2" == false && "$RUN_LASSOSUM2" == false ]]; then
    usage
fi



# --- METHOD 1: Clumping + Thresholding (C+T) ---
if [[ "$RUN_CT" == true ]]; then
    echo "[$(date)] Starting C+T Pipeline..."
    mkdir -p ${output_path}/gwas
    
    # 1. File Generation
    Rscript "${path_repo}/src/prepare_sumstats.R" \
        --input  ${summary_stats_file} \
        --bim "$bim_file_path" \
        --n_total "$n_total_gwas" \
        --output "${output_path}/gwas/CT_PRSice2_summary_stat_file.txt"
    
    # 2. Config & Run
    mkdir -p ${output_path}/prs_pipeline/CT/temp
    CT_CONFIG="${output_path}/prs_pipeline/CT/temp/CT_temp_config.txt"
    awk '{print $1, $2, $6}' OFS="\t" ${study_sample}.fam > ${output_path}/gwas/study_sample_pheno.txt
    cat <<EOF > "$CT_CONFIG"
study_sample=${study_sample}
sum_stats_file=${output_path}/gwas/CT_PRSice2_summary_stat_file.txt
phenotype_info_file=${output_path}/gwas/study_sample_pheno.txt
gwas_pca_eigenvec_file=${gwas_pca_eigenvec_file}
output_path=${output_path}
path_prs_pipeline=${path_repo}
EOF

    bash "${path_repo}/src/run_CT.sh" --c "$CT_CONFIG"
fi

# --- METHOD 2: LDpred2 ---
if [[ "$RUN_LDPRED2" == true ]]; then
    echo "[$(date)] Starting LDpred2 Pipeline..."
    mkdir -p "${output_path}/prs_pipeline/LDpred2"
    echo "Running below
    Rscript ${path_repo}/src/run_LDpred2.R \
        --anc_bed ${study_sample}.bed \
        --ss $qassoc_path \
        --bim $bim_file_path \
        --beta_se $table_with_beta_se_path \
        --out ${output_path}/prs_pipeline/LDpred2/prs_method
        "
    Rscript "${path_repo}/src/run_LDpred2.R" \
        --anc_bed "${study_sample}.bed" \
        --ss "$qassoc_path" \
        --bim "$bim_file_path" \
        --beta_se "$table_with_beta_se_path" \
        --out "${output_path}/prs_pipeline/LDpred2/prs_method"
fi

# --- METHOD 3: lassosum2 ---
if [[ "$RUN_LASSOSUM2" == true ]]; then
    echo "[$(date)] Starting lassosum2 Pipeline..."
    mkdir -p "${output_path}/prs_pipeline/lassosum2"
    
    echo "Running below
    Rscript ${path_repo}/src/run_lassosrum2.R \
        --anc_bed ${study_sample}.bed \
        --ss $qassoc_path \
        --bim $bim_file_path \
        --beta_se $table_with_beta_se_path \
        --out ${output_path}/prs_pipeline/lassosum2/prs_method
        "
    # Note: Fixed the typo from your original script (calling run_LDpred2.R instead of run_lassosum2.R)
    Rscript "${path_repo}/src/run_lassosum2.R" \
        --anc_bed "${study_sample}.bed" \
        --ss "$qassoc_path" \
        --bim "$bim_file_path" \
        --beta_se "$table_with_beta_se_path" \
        --out "${output_path}/prs_pipeline/lassosum2/prs_method"
fi

echo "All requested PRS methods have completed."
