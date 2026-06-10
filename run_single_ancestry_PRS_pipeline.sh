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
gwas_pca_eigenvec_file=""
afreq_file=""
ncores=16
ld_cache_dir=""
ld_matrix_dir=""
skip_ss_generation=0
binary_flag=F # accepts T/F
phenotype_info_file=""
test_sample=""
test_pca_eigenvec_file=""

# Multi-phenotype parameters
summary_stats_files=""       # Comma-separated list of GWAS files, one per phenotype
multi_pheno_file=""          # Tab-sep file: FID, IID, pheno1, pheno2, ... (with header)

# --- Environment ---
# (Environment loaded via Singularity container - conda activation removed)

# --- Pre-load config ---
for arg in "$@"; do
    if [[ "$arg" == "-C" ]]; then
        # Extract value after -C (macOS-compatible, no -P flag)
        conf_file=$(echo "$@" | sed -n 's/.*-C \([^ ]*\).*/\1/p')
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
        C) ;; # Handled in pre-load step
        c) RUN_CT=true ;;
        l) RUN_LDPRED2=true ;;
        s) RUN_LASSOSUM2=true ;;
        P) RUN_PRSice2=true ;;
        S) skip_ss_generation=1 ;;
        B) binary_flag=T ;;
        *) usage ;;
    esac
done

if [[ "$RUN_CT" == false && "$RUN_LDPRED2" == false && "$RUN_LASSOSUM2" == false && "$RUN_PRSice2" == false ]]; then
    usage
fi

# ===================================================================
# PIPELINE BODY (shared by single- and multi-phenotype modes)
# ===================================================================
# Usage: run_phenotype_pipeline <ss_file> <pheno_name>
#   ss_file    — GWAS summary stats for this phenotype
#   pheno_name — phenotype name (used for output subdir & column lookup);
#                empty string = single-phenotype mode (backward compat)
# ===================================================================
run_phenotype_pipeline() {
    local ss_file="$1"
    local pheno_name="$2"

    local tag=""
    if [[ -n "$pheno_name" ]]; then
        tag="[${pheno_name}]"
    fi

    # Output base for this phenotype's method directories
    local methods_output="${output_path}/prs_pipeline"
    if [[ -n "$pheno_name" ]]; then
        methods_output="${methods_output}/${pheno_name}"
    fi

    echo "[$(date)] $tag Starting pipeline (methods_output=${methods_output})"

    # --- 1. Summary Stats Prep ---
    local ss_local="$ss_file"
    mkdir -p "${output_path}/gwas"
    if [[ "$skip_ss_generation" == 0 ]]; then
        echo "[$(date)] $tag Aligning summary stats..."
        Rscript "${path_repo}/src/prepare_sumstats.R" \
            --input "$ss_local" \
            --bim "$bim_file_path" \
            --n_total "$n_total_gwas" \
            --output "${output_path}/gwas/CT_PRSice2_summary_stat_file.txt"
        ss_local="${output_path}/gwas/CT_PRSice2_summary_stat_file.txt"
    fi

    # --- 2. Phenotype File ---
    if [[ -n "$multi_pheno_file" && -f "$multi_pheno_file" && -n "$pheno_name" ]]; then
        # Multi-pheno mode: extract the named column
        echo "[$(date)] $tag Extracting phenotype '${pheno_name}' from multi-pheno file"
        # Find column index by header name
        local pheno_col
        pheno_col=$(head -1 "$multi_pheno_file" | tr '\t' '\n' | awk -v name="$pheno_name" 'NR==1{next} $0==name{print NR}')
        if [[ -z "$pheno_col" ]]; then
            echo "Error: Phenotype column '${pheno_name}' not found in header of ${multi_pheno_file}" >&2
            echo "  Header columns: $(head -1 "$multi_pheno_file" | tr '\t' ',')" >&2
            exit 1
        fi
        echo "[$(date)] $tag   Column ${pheno_col}: ${pheno_name}"
        awk -v col="$pheno_col" 'NR==1{print $1, $2, $col} NR>1{print $1, $2, $col}' OFS="\t" "$multi_pheno_file" \
            > "${output_path}/gwas/study_sample_pheno.txt"
    elif [[ -n "${phenotype_info_file:-}" && -f "$phenotype_info_file" ]]; then
        echo "[$(date)] $tag Using external phenotype file: $phenotype_info_file"
        # Auto-detect and add header if missing (check first field of line 1)
        pheno_header=$(head -1 "$phenotype_info_file" | awk '{print $1}')
        if [[ "$pheno_header" != "FID" && "$pheno_header" != "fid" ]]; then
            ncols=$(head -1 "$phenotype_info_file" | awk '{print NF}')
            echo "[$(date)] $tag   No header detected — prepending ${ncols}-column header (FID, IID, phenotype1...phenotype$((ncols-2)))"
            {
                awk -v ncols="$ncols" 'BEGIN{
                    printf "FID\tIID"
                    for(i=3; i<=ncols; i++) printf "\tphenotype%d", i-2
                    print ""
                }'
                cat "$phenotype_info_file"
            } > "${output_path}/gwas/study_sample_pheno.txt"
        else
            cp "$phenotype_info_file" "${output_path}/gwas/study_sample_pheno.txt"
        fi
    else
        echo "[$(date)] $tag Extracting phenotype from .fam file (column 6) with header"
        awk 'BEGIN{print "FID\tIID\tphenotype"} {print $1, $2, $6}' OFS="\t" "${study_sample}.fam" \
            > "${output_path}/gwas/study_sample_pheno.txt"
    fi

    # --- 3. LD Matrix ---
    local ld_dir_local="$ld_matrix_dir"
    if [[ -n "$ld_matrix_dir" && ( "$RUN_LDPRED2" == true || "$RUN_LASSOSUM2" == true ) ]]; then
        if [[ -n "$pheno_name" ]]; then
            ld_dir_local="${ld_matrix_dir}/${pheno_name}"
        fi
        mkdir -p "$ld_dir_local"
        if [[ ! -f "${ld_dir_local}/map.rds" ]]; then
            echo "[$(date)] $tag Generating LD matrix..."
            Rscript "${path_repo}/src/generate_ld_matrix.R" \
                --anc_bed "${study_sample}.bed" \
                --out "$ld_dir_local" \
                --ncores "$ncores"
        else
            echo "[$(date)] $tag Using existing LD matrix"
        fi
    fi

    # --- 4. C+T ---
    if [[ "$RUN_CT" == true ]]; then
        (
        echo "[$(date)] $tag Starting C+T Pipeline..."
        mkdir -p ${methods_output}/CT/temp
        CT_CONFIG="${methods_output}/CT/temp/CT_temp_config.txt"
        cat <<EOF > "$CT_CONFIG"
study_sample=${study_sample}
sum_stats_file=${ss_local}
phenotype_info_file=${output_path}/gwas/study_sample_pheno.txt
output_path=${methods_output}
path_prs_pipeline=${path_repo}
EOF
        # Only pass PCA file if user explicitly provided one
        if [[ -n "${gwas_pca_eigenvec_file:-}" ]]; then
            echo "gwas_pca_eigenvec_file=${gwas_pca_eigenvec_file}" >> "$CT_CONFIG"
        fi
        bash "${path_repo}/src/run_CT.sh" --c "$CT_CONFIG"
        ) &
    fi

    # --- 5. LDpred2 ---
    if [[ "$RUN_LDPRED2" == true ]]; then
        (
        echo "[$(date)] $tag Starting LDpred2 Pipeline..."
        mkdir -p "${methods_output}/LDpred2"
        LDpred2_args="--anc_bed ${study_sample}.bed --ss $ss_local --bim $bim_file_path --out ${methods_output}/LDpred2/prs_method --ncores $ncores --pheno ${output_path}/gwas/study_sample_pheno.txt"
        if [[ -n "$afreq_file" ]]; then
            LDpred2_args="$LDpred2_args --afreq $afreq_file"
        fi
        if [[ -n "$ld_cache_dir" ]]; then
            LDpred2_args="$LDpred2_args --ld-cache-dir $ld_cache_dir"
        fi
        if [[ -n "$ld_dir_local" ]]; then
            LDpred2_args="$LDpred2_args --ld-matrix-dir $ld_dir_local"
        fi
        Rscript "${path_repo}/src/run_LDpred2.R" $LDpred2_args
        ) &
    fi

    # --- 6. lassosum2 ---
    if [[ "$RUN_LASSOSUM2" == true ]]; then
        (
        echo "[$(date)] $tag Starting lassosum2 Pipeline..."
        mkdir -p "${methods_output}/lassosum2"
        lassosum2_args="--anc_bed ${study_sample}.bed --ss $ss_local --bim $bim_file_path --out ${methods_output}/lassosum2/prs_method --ncores $ncores --pheno ${output_path}/gwas/study_sample_pheno.txt"
        if [[ -n "$afreq_file" ]]; then
            lassosum2_args="$lassosum2_args --afreq $afreq_file"
        fi
        if [[ -n "$ld_cache_dir" ]]; then
            lassosum2_args="$lassosum2_args --ld-cache-dir $ld_cache_dir"
        fi
        if [[ -n "$ld_dir_local" ]]; then
            lassosum2_args="$lassosum2_args --ld-matrix-dir $ld_dir_local"
        fi
        Rscript "${path_repo}/src/run_lassosum2.R" $lassosum2_args
        ) &
    fi

    # --- 7. PRSice2 ---
    if [[ "$RUN_PRSice2" == true ]]; then
        echo "[$(date)] $tag Starting PRSice2 Pipeline..."
        local prsice_out="${methods_output}/PRSice2/prs_method"
        mkdir -p "$prsice_out"
        bash "${path_repo}/src/run_PRSice2.sh" \
            "$ss_local" \
            "$study_sample" \
            "$binary_flag" \
            "${output_path}/gwas/study_sample_pheno.txt" \
            "${methods_output}" \
            "${path_repo}" \
            "${prsice_out}"
    fi

    echo "[$(date)] $tag All jobs submitted to background. Waiting..."
    wait

    # --- 8. Test Evaluation (if test data provided) ---
    if [[ -n "${test_sample:-}" ]]; then
        echo "[$(date)] $tag Starting test evaluation..."
        bash "${path_repo}/src/score_test.sh" \
            --test-bfile "$test_sample" \
            --train-out-dir "$methods_output" \
            --pheno-file "${phenotype_info_file:-}" \
            --test-pca-file "${test_pca_eigenvec_file:-}" \
            --sumstats "$ss_local" \
            --path-repo "$path_repo" \
            --binary-flag "$binary_flag" \
            --ran-ct "$RUN_CT" \
            --ran-ldpred2 "$RUN_LDPRED2" \
            --ran-lassosum2 "$RUN_LASSOSUM2" \
            --ran-prsice2 "$RUN_PRSice2"
        echo "[$(date)] $tag Test evaluation complete."
    fi

    echo "[$(date)] $tag Pipeline complete."
}

# ===================================================================
# DISPATCH: single-phenotype vs multi-phenotype
# ===================================================================

# Determine whether we are in multi-phenotype mode.
# Multi-pheno is active when summary_stats_files is non-empty AND
# either contains a comma (multiple files) OR multi_pheno_file is set.
multi_mode=false
if [[ -n "$summary_stats_files" ]]; then
    # Count comma-separated entries
    IFS=',' read -ra SS_FILES <<< "$summary_stats_files"
    if [[ ${#SS_FILES[@]} -ge 2 ]]; then
        multi_mode=true
    elif [[ ${#SS_FILES[@]} -eq 1 && -n "$multi_pheno_file" ]]; then
        multi_mode=true
        # Single file + multi_pheno_file still triggers multi-mode
        # (user might want named output subdirs even for 1 pheno)
    fi
fi

if [[ "$multi_mode" == true ]]; then
    echo "[$(date)] Multi-phenotype mode enabled"
    echo "  summary_stats_files: ${summary_stats_files}"
    echo "  multi_pheno_file: ${multi_pheno_file}"

    # Validate multi_pheno_file
    if [[ -z "$multi_pheno_file" || ! -f "$multi_pheno_file" ]]; then
        echo "Error: multi_pheno_file is required in multi-phenotype mode" >&2
        exit 1
    fi

    # Read phenotype names from the multi_pheno_file header (columns 3+)
    # Ensure tab-delimited columns are parsed correctly
    pheno_names=()
    IFS=$'\t' read -ra header_cols < <(head -1 "$multi_pheno_file")
    for ((i=2; i<${#header_cols[@]}; i++)); do
        pheno_names+=("${header_cols[$i]}")
    done

    echo "  Phenotypes from header: ${pheno_names[*]}"
    echo "  Number of summary stats files: ${#SS_FILES[@]}"
    echo "  Number of phenotype columns: ${#pheno_names[@]}"

    if [[ ${#SS_FILES[@]} -ne ${#pheno_names[@]} ]]; then
        echo "Error: Number of summary stats files (${#SS_FILES[@]}) does not match" >&2
        echo "  number of phenotype columns (${#pheno_names[@]})" >&2
        echo "  Summary stats files: ${SS_FILES[*]}" >&2
        echo "  Phenotype columns: ${pheno_names[*]}" >&2
        exit 1
    fi

    # Run each phenotype sequentially (methods within each run in parallel)
    for ((i=0; i<${#SS_FILES[@]}; i++)); do
        ss="${SS_FILES[$i]}"
        pname="${pheno_names[$i]}"
        echo ""
        echo "============================================================"
        echo "[$(date)] Processing phenotype ${i}: ${pname}"
        echo "  Summary stats: ${ss}"
        echo "============================================================"
        run_phenotype_pipeline "$ss" "$pname"
        echo "[$(date)] Finished phenotype ${i}: ${pname}"
        echo ""
    done

    echo "[$(date)] All phenotypes completed."
else
    # --- Original single-phenotype behavior ---
    # Use summary_stats_file (singular) as fallback if summary_stats_files not set
    local_ss="${summary_stats_file}"
    if [[ -n "$summary_stats_files" ]]; then
        # Single entry in summary_stats_files; parse it
        IFS=',' read -ra SS_FILES <<< "$summary_stats_files"
        local_ss="${SS_FILES[0]}"
    fi
    run_phenotype_pipeline "$local_ss" ""
fi
