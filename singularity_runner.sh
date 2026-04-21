#!/bin/bash
# Wrapper script to run PRS pipeline with Singularity container
# Works on both SLURM HPC and local machines

set -eu

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF_PATH="${SCRIPT_DIR}/singleprs_latest.sif"
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
        --bind "${SCRIPT_DIR}:${SCRIPT_DIR}" \
        --bind "${binds}" \
        --bind "/tmp" \
        --pwd "${SCRIPT_DIR}" \
        "${SIF_PATH}" \
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
                    local mount_point="${path%%/*}"
                    if [[ -n "$mount_point" && ! "$binds" =~ "$mount_point" ]]; then
                        binds="${binds:+${binds},}${mount_point}"
                    fi
                    ;;
            esac
        done < "${CONFIG_FILE}"
    fi
    
    echo "$binds"
}

# Parse own args first to find config
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bind) BIND_PATHS="$2"; shift 2 ;;
        --) shift; break ;;
        --C|-C) CONFIG_FILE="$2"; shift 2 ;;
        -*) shift ;;
        *) break ;;
    esac
done

# Collect any remaining args for pipeline
PIPELINE_ARGS="$@"

# Build bind paths
BIND_PATHS="${BIND_PATHS:-$(auto_bind)}"
if [[ -n "${CONFIG_FILE:-}" && ! "$BIND_PATHS" =~ "/projects" ]]; then
    BIND_PATHS="${BIND_PATHS:+${BIND_PATHS},}/projects"
fi

# Execute in container
run_in_container "${BIND_PATHS}" \
    bash "${SCRIPT_DIR}/run_single_ancestry_PRS_pipeline.sh" ${PIPELINE_ARGS}