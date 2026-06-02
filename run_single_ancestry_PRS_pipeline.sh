#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH -p agsmall
#SBATCH -o prs_unified_%j.out
#SBATCH --job-name prs_pipeline

set -eu

# --- Defaults ---
RUN_CT=false
RUN_LDPRED2=false
RUN_LASSOSUM2=false
RUN_PRSice2=false

# Paths (Consolidated)
summary_stats_file="/projects/standard/gdc/public/prs_methods/data/Phenotypes/updated_summary_stats/updated_phenocode-aseg_lh_volume_Hippocampus.tsv"
bim_file_path="/projects/standard/gdc/public/prs_methods/data/simulated_1000G/anc1_plink_files/archived/AFR_simulation_gwas.bim"
study_sample="/projects/standard/gdc/public/prs_methods/data/simulated_1000G/anc1_plink_files/AFR_simulation_study_sample"
output_path="/projects/standard/gdc/public/prs_methods/data/simulated_1000G"
path_repo="/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline"
n_total_gwas=31968
gwas_pca_eigenvec_file="/projects/standard/gdc/shared/abcdTest/04-globalAncestry/merged_dataset_pca.eigenvec"
afreq_file=""
ncores=16
ld_cache_dir=""
ld_matrix_dir=""
skip_ss_generation=0
binary_flag=F # accepts T/F
phenotype_info_file=""

# --- Environment ---
# (Environment loaded via Singularity container - conda activation removed)

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
    echo "  -P    Run PRSice2. *For this to run as intended (C+T) needs to be run as well.*"
    echo "  -S    Skip summary stats allignment step. Not recommended."
    echo "  -B    Include if you are using a binary phenotype."
    exit 1
}

# --- Arg Parser ---
while getopts "C:clsPSB" opt; do
    case "$opt" in
        C) ;; # Handled in pre-load step, but kept here so getopts doesn't complain
        c) RUN_CT=true ;;
        l) RUN_LDPRED2=true ;;
        s) RUN_LASSOSUM2=true ;;
        P) RUN_PRSice2=true ;;
        S) skip_ss_generation=1 ;;
        B) binary_flag=T ;;
        *) usage ;;
    esac
done

# If no flags provided, show usage
if [[ "$RUN_CT" == false && "$RUN_LDPRED2" == false && "$RUN_LASSOSUM2" == false && "$RUN_PRSice2" == false ]]; then
    usage
fi

if [[ "$skip_ss_generation" == 0 ]]; then
  mkdir -p ${output_path}/gwas
    
  # 1. File Generation
  Rscript "${path_repo}/src/prepare_sumstats.R" \
    --input  ${summary_stats_file} \
    --bim "$bim_file_path" \
    --n_total "$n_total_gwas" \
    --output "${output_path}/gwas/CT_PRSice2_summary_stat_file.txt"
  summary_stats_file="${output_path}/gwas/CT_PRSice2_summary_stat_file.txt"
fi

if [[ -n "${phenotype_info_file:-}" && -f "$phenotype_info_file" ]]; then
  echo "[$(date)] Using external phenotype file: $phenotype_info_file"
  cp "$phenotype_info_file" "${output_path}/gwas/study_sample_pheno.txt"
else
  echo "[$(date)] Extracting phenotype from .fam file (column 6)"
  awk '{print $1, $2, $6}' OFS="\t" ${study_sample}.fam > ${output_path}/gwas/study_sample_pheno.txt
fi

# --- LD MATRIX GENERATION (prep step for LDpred2 / lassosum2) ---
if [[ -n "$ld_matrix_dir" && ( "$RUN_LDPRED2" == true || "$RUN_LASSOSUM2" == true ) ]]; then
  mkdir -p "$ld_matrix_dir"
  if [[ ! -f "${ld_matrix_dir}/map.rds" ]]; then
    echo "[$(date)] Generating pre-computed LD matrix..."
    Rscript "${path_repo}/src/generate_ld_matrix.R" \
      --anc_bed "${study_sample}.bed" \
      --out "$ld_matrix_dir" \
      --ncores "$ncores"
    echo "[$(date)] LD matrix saved to $ld_matrix_dir"
  else
    echo "[$(date)] Using existing LD matrix from $ld_matrix_dir"
  fi
fi

# --- METHOD 1: Clumping + Thresholding (C+T) ---
if [[ "$RUN_CT" == true ]]; then
   ( 
    echo "[$(date)] Starting C+T Pipeline..."
    
    # 2. Config & Run
  mkdir -p ${output_path}/prs_pipeline/CT/temp
  CT_CONFIG="${output_path}/prs_pipeline/CT/temp/CT_temp_config.txt"
  cat <<EOF > "$CT_CONFIG"
study_sample=${study_sample}
sum_stats_file=${summary_stats_file}
phenotype_info_file=${output_path}/gwas/study_sample_pheno.txt
gwas_pca_eigenvec_file=${gwas_pca_eigenvec_file}
output_path=${output_path}
path_prs_pipeline=${path_repo}
EOF

  bash "${path_repo}/src/run_CT.sh" --c "$CT_CONFIG"
  ) &   
fi

# --- METHOD 2: LDpred2 ---
if [[ "$RUN_LDPRED2" == true ]]; then
    (
    echo "[$(date)] Starting LDpred2 Pipeline..."
    mkdir -p "${output_path}/prs_pipeline/LDpred2"
    echo "[DEBUG] afreq_file='${afreq_file:-}'" >&2
    LDpred2_args="--anc_bed ${study_sample}.bed --ss $summary_stats_file --bim $bim_file_path --out ${output_path}/prs_pipeline/LDpred2/prs_method --ncores $ncores --pheno ${output_path}/gwas/study_sample_pheno.txt"
    if [[ -n "$afreq_file" ]]; then
      LDpred2_args="$LDpred2_args --afreq $afreq_file"
    fi
    if [[ -n "$ld_cache_dir" ]]; then
      LDpred2_args="$LDpred2_args --ld-cache-dir $ld_cache_dir"
    fi
    if [[ -n "$ld_matrix_dir" ]]; then
      LDpred2_args="$LDpred2_args --ld-matrix-dir $ld_matrix_dir"
    fi
    echo "Running below
    Rscript ${path_repo}/src/run_LDpred2.R $LDpred2_args"
    Rscript "${path_repo}/src/run_LDpred2.R" $LDpred2_args
    ) &
fi

# --- METHOD 3: lassosum2 ---
if [[ "$RUN_LASSOSUM2" == true ]]; then
    (
    echo "[$(date)] Starting lassosum2 Pipeline..."
    mkdir -p "${output_path}/prs_pipeline/lassosum2"
    echo "[DEBUG] afreq_file='${afreq_file:-}'" >&2
    
    lassosum2_args="--anc_bed ${study_sample}.bed --ss $summary_stats_file --bim $bim_file_path --out ${output_path}/prs_pipeline/lassosum2/prs_method --ncores $ncores"
    if [[ -n "$afreq_file" ]]; then
      lassosum2_args="$lassosum2_args --afreq $afreq_file"
    fi
    if [[ -n "$ld_cache_dir" ]]; then
      lassosum2_args="$lassosum2_args --ld-cache-dir $ld_cache_dir"
    fi
    if [[ -n "$ld_matrix_dir" ]]; then
      lassosum2_args="$lassosum2_args --ld-matrix-dir $ld_matrix_dir"
    fi
    echo "Running below
    Rscript ${path_repo}/src/run_lassosum2.R $lassosum2_args"

    Rscript "${path_repo}/src/run_lassosum2.R" $lassosum2_args
    ) &
fi


# --- METHOD 4: PRSice2 ---
if [[ "$RUN_PRSice2" == true ]]; then
    
    echo "[$(date)] Starting PRSice2 Pipeline..."
    mkdir -p "${output_path}/prs_pipeline/PRSice2"
    
    echo "Running below
    bash ${path_repo}/src/run_PRSice2.sh \
        $summary_stats_file \
        $study_sample \
        $binary_flag \
        ${output_path}/gwas/study_sample_pheno.txt \
        ${output_path} \
        ${path_repo}
        "

# phenotype_file has FID IID phenotype # as a column header 
    bash "${path_repo}/src/run_PRSice2.sh" \
        "$summary_stats_file" \
        "$study_sample" \
        "$binary_flag" \
        "${output_path}/gwas/study_sample_pheno.txt" \
        "${output_path}" \
        "${path_repo}"
    
fi

echo "All jobs submitted to background. Waiting for completion of all chosen methods."
wait
echo "All requested PRS methods have completed."
