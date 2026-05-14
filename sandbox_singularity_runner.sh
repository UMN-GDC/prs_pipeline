#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH -p small
#SBATCH -o prs_unified_%j.out
#SBATCH --job-name prs_pipeline

# Wrapper script to run PRS pipeline with Singularity container
# Works on both SLURM HPC and local machines

set -eu

#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#SIF_PATH="${SCRIPT_DIR}/prsv2_latest.sif"
ENV_NAME="singlePRS"

usage() {
    echo "Usage: $0 [singularity options] [--C config_file] [pipeline options]"
    echo ""
    echo "Singularity options:"
    echo "  --bind <paths>    Bind paths (default: auto-detect from config)"
    echo ""
    echo "Pipeline options:"
    echo "  -C <config>     Path to config file"
    echo "  -c             Run Clumping + Thresholding (C+T)"
    echo "  -l             Run LDpred2"
    echo "  -s             Run lassosum2"
    echo "  -P             Run PRSice2 (requires -c)"
    echo "  -S             Skip summary stats alignment"
    echo "  -B             Binary phenotype"
    echo ""
    echo "Examples:"
    echo "  # On SLURM HPC"
    echo "  sbatch $0 --C config.txt -c -P"
    echo ""
    echo "  # Local machine"
    echo "  $0 --C config.txt -c -l -s"
    exit 1
}

is_slurm() {
    [[ -n "${SLURM_JOB_ID:-}" ]]
}

run_in_container() {
    local binds="$1"
    shift
    
    singularity exec \
        ${binds:+--bind "$binds"} \
        --bind "/tmp:/tmp" \
        --pwd "${path_repo}" \
        "${path_repo}/prsv2_latest.sif" \
        "$@"
}

auto_bind() {
    local binds=""
    
    if [[ -f "${CONFIG_FILE:-}" ]]; then
        while IFS= read -r line; do
            case "$line" in
                summary_stats_file=*|bim_file_path=*|study_sample=*|output_path=*|path_repo=*|gwas_pca_eigenvec_file=*)
                    local path="${line#*=}"
                    path="${path%\"}"
                    path="${path#\"}"
                    
                    if [[ "$path" == /* ]]; then
                        # Get the top-level directory (e.g., /scratch.global)
                        local root_dir="/$(echo "$path" | cut -d'/' -f2)"
                        
                        # CRITICAL: Only add if the directory exists on the host MSI machine
                        if [[ -d "$root_dir" ]]; then
                            if [[ ! "$binds" =~ "$root_dir" ]]; then
                                binds="${binds:+${binds},}${root_dir}"
                            fi
                        else
                            # Optional: warn the user if their config points to a ghost path
                            echo "Note: Config references $root_dir, but it doesn't exist on this system. Skipping bind." >&2
                        fi
                    fi
                    ;;
            esac
        done < "${CONFIG_FILE}"
    fi
    echo "$binds"
}

PIPELINE_ARGS=""
BIND_PATHS=""

# Parse own args first to find config
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bind) 
            BIND_PATHS="$2"; shift 2 ;;
        --C|-C) 
            CONFIG_FILE=$(realpath "$2"); 
            source ${CONFIG_FILE}
            PIPELINE_ARGS+="-C $2 "
            shift 2 ;;
        --) 
            shift; PIPELINE_ARGS+="$* "; break ;;
        -*) 
            PIPELINE_ARGS+="$1 "; shift ;;
        *) 
            PIPELINE_ARGS+="$1"; shift ;;
    esac
done


# Build bind paths only from existing directories
DETECTED_BINDS="$(auto_bind)"
FINAL_BINDS="${BIND_PATHS:-$DETECTED_BINDS}"

# Execute
run_in_container "${FINAL_BINDS}" \
    bash "${path_repo}/run_single_ancestry_PRS_pipeline.sh" ${PIPELINE_ARGS}
