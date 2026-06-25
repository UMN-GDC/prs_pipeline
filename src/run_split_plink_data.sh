#!/bin/bash 
#SBATCH --time=1:00:00
#SBATCH --ntasks=8
#SBATCH --mem=16g
#SBATCH --tmp=32g
#SBATCH -p agsmall
#SBATCH -o split_plink_data.out
#SBATCH --job-name split_plink_data

# set -euo pipefail
# shopt -s extglob

### --- CONFIGURABLE DEFAULTS -------------------------------------------------
# Defaults (can be overridden via CLI)
plink_file_anc1="/projects/standard/gdc/public/prs_methods/data/test/sim_1/AFR_simulation_study_sample"                     # target ancestry
plink_file_anc2="/projects/standard/gdc/public/prs_methods/data/test/sim_1/EUR_simulation_study_sample"                     # training ancestry

train_percent=50
valid_percent=0
test_percent=50
rand_seed=42
no_plink=0
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located
# path_R_packages=/projects/standard/gdc/public/Ref/R
# export R_LIBS_USER="$path_R_packages"


### --- FUNCTIONS ------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Split PLINK files into train/val/test subsets. When validation percent is 0
(default), produces a 2-way train/test split instead.

Options:
  -1 <plink_file_anc1>      Full path to target ancestry plink files (default: ${plink_file_anc1})
  -2 <plink_file_anc2>      Full path to training ancestry plink files  (default: ${plink_file_anc2})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -t <train_percent>        Percent of data for training (default: 50)
  -v <validation_percent>   Percent of data for validation (default: 0; set >0 for 3-way split)
  -T <test_percent>         Percent of data for testing (default: 50)
  -S <seed>                 Randomization seed (default: 42)
  -N                        Skip generating PLINK separated files (sample lists only)
  -h                        Show this help and exit

Examples:
  2-way split (train/test, default):
    bash $0 -1 /path/to/AFR_study_sample -2 /path/to/EUR_study_sample

  3-way split (train/val/test):
    bash $0 -1 /path/to/AFR_study_sample -2 /path/to/EUR_study_sample -t 70 -v 15 -T 15

  Sample lists only (no PLINK file generation):
    bash $0 -1 /path/to/AFR_study_sample -2 /path/to/EUR_study_sample -N
EOF
  exit 1
}

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}



### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:r:t:T:v:S:Nh" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    t) train_percent="$OPTARG" ;;
    T) test_percent="$OPTARG" ;;
    v) valid_percent="$OPTARG" ;;
    S) rand_seed="$OPTARG" ;;
    N) no_plink=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done


# Load environments / modules
# module load R/4.4.0-openblas-rocky8
module load plink
module load plink/2.00-alpha-091019

base_location=$(dirname "$plink_file_anc1")
mkdir -p "${base_location}/randomization_ids_anc1"
mkdir -p "${base_location}/randomization_ids_anc2"

# Split phenotype / plink files
source /projects/standard/gdc/public/envs/load_miniconda3.sh
log "Splitting plink data for "$plink_file_anc1""
if [ "${no_plink}" -eq 1 ]; then
  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed} \
    --no_plink

  mv "${base_location}"/*samples.txt "${base_location}/randomization_ids_anc1"
  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed} \
    --no_plink
  mv "${base_location}"/*samples.txt "${base_location}/randomization_ids_anc2"
else

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}
  mv "${base_location}"/*samples.txt "${base_location}/randomization_ids_anc1"

  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}
  mv "${base_location}"/*samples.txt "${base_location}/randomization_ids_anc2"
fi

conda deactivate
conda deactivate


if [[ "${valid_percent}" -gt 0 ]]; then
  echo "Finished generating inputs (3-way split): train=${train_percent}%, val=${valid_percent}%, test=${test_percent}%"
else
  echo "Finished generating inputs (2-way split): train=${train_percent}%, test=${test_percent}%"
fi

