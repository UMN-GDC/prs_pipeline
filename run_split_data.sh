#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH -p medium
#SBATCH -o split_data_%j.out
#SBATCH --job-name split_data

# Wrapper to split PLINK data into train/val/test subsets inside the
# prsv2_latest.sif container. All dependencies (plink, plink2, python)
# live inside the container — no host modules needed.
#
# Accepts one or two PLINK inputs and forwards them to
# src/run_split_plink_data.sh after launching the container.
#
# Examples:
#   sbatch run_split_data.sh -1 /path/to/study_sample -r /path/to/prs_pipeline
#   bash run_split_data.sh -1 /path/to/anc1 -2 /path/to/anc2 -r /path/to/prs_pipeline -t 70 -T 30

set -eu

usage() {
    cat <<EOF
Usage: $0 [options]

Split PLINK data into train/test (or train/val/test) subsets inside the
Singularity container. All PLINK and Python dependencies are provided by
prsv2_latest.sif — no host modules needed.

Options: (see src/run_split_plink_data.sh for full details)
  -1 <path>    PLINK prefix (required)
  -2 <path>    Second PLINK prefix (optional)
  -r <path>    Path to prs_pipeline repository (required)
  -t <pct>     Training percent (default: 50)
  -v <pct>     Validation percent (default: 0; 0 = 2-way train/test)
  -T <pct>     Test percent (default: 50)
  -S <int>     Random seed (default: 42)
  -N           Sample lists only (skip PLINK file generation)
  -h           Show this help

Examples:
  Single ancestry (1 input):
    sbatch $0 -1 /path/to/study_sample -r /path/to/prs_pipeline

  Two ancestries (2 inputs):
    bash $0 -1 /path/to/AFR -2 /path/to/EUR -r /path/to/prs_pipeline -t 70 -T 30 -N
EOF
    exit 1
}

# Parse args (same flags as run_split_plink_data.sh)
train_percent=50
valid_percent=0
test_percent=50
rand_seed=42
no_plink=0
plink_file_anc1=""
plink_file_anc2=""
path_to_repo=""

while getopts ":1:2:r:t:v:T:S:Nh" opt; do
    case "$opt" in
        1) plink_file_anc1="$OPTARG" ;;
        2) plink_file_anc2="$OPTARG" ;;
        r) path_to_repo="$OPTARG" ;;
        t) train_percent="$OPTARG" ;;
        v) valid_percent="$OPTARG" ;;
        T) test_percent="$OPTARG" ;;
        S) rand_seed="$OPTARG" ;;
        N) no_plink=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$plink_file_anc1" || -z "$path_to_repo" ]]; then
    echo "ERROR: -1 and -r are required"
    usage
fi

SIF="${path_to_repo}/prsv2_latest.sif"
if [[ ! -f "$SIF" ]]; then
    echo "ERROR: Container not found at $SIF"
    echo "  Pull it with: apptainer pull oras://ghcr.io/mainsqu33ze/gdcgenomicsqc/prsv2:latest"
    exit 1
fi

# Auto-detect bind mounts: extract the top-level root from each path
# and add it if the directory exists on the host.
auto_bind() {
    local binds=""
    for path in "$plink_file_anc1" "$path_to_repo" ${plink_file_anc2:+"$plink_file_anc2"}; do
        local dir
        dir=$(dirname "$path")
        local root
        root="/$(echo "$dir" | cut -d'/' -f2)"
        if [[ -d "$root" && ! "$binds" =~ "$root" ]]; then
            binds="${binds:+${binds},}${root}"
        fi
    done
    echo "$binds"
}

BINDS=$(auto_bind)

echo "[$(date)] Splitting data inside container: $(basename "$SIF")"
echo "  Binds: $BINDS"
echo "  Split: train=${train_percent}% val=${valid_percent}% test=${test_percent}%"

# Build forwarded args — only include -2 when provided
FORWARD_ARGS="-1 $plink_file_anc1"
if [[ -n "$plink_file_anc2" ]]; then
    FORWARD_ARGS="$FORWARD_ARGS -2 $plink_file_anc2"
fi
FORWARD_ARGS="$FORWARD_ARGS -r $path_to_repo -t $train_percent -v $valid_percent -T $test_percent -S $rand_seed"
if [[ "$no_plink" -eq 1 ]]; then
    FORWARD_ARGS="$FORWARD_ARGS -N"
fi

singularity exec \
    --bind "$BINDS" \
    --bind "/tmp:/tmp" \
    --pwd "${path_to_repo}" \
    "$SIF" \
    bash "${path_to_repo}/src/run_split_plink_data.sh" \
        $FORWARD_ARGS
