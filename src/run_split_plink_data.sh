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

gwas_percent=0
train_percent=50
valid_percent=20
test_percent=30
rand_seed=42
no_plink=0
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located
# path_R_packages=/projects/standard/gdc/public/Ref/R
# export R_LIBS_USER="$path_R_packages"


### --- FUNCTIONS ------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <plink_file_anc1>      Full path to target ancestry plink files (default: ${plink_file_anc1})
  -2 <plink_file_anc2>      Full path to training ancestry plink files  (default: ${plink_file_anc2})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -t <train_percent>        Percent of data for training (default: 50)
  -T <testing_percent>      Percent of data for testing (default: 20)
  -v <validation_percent>   Percent of data for validation (default: 30)
  -S <seed>                 Randomization seed (default:42)
  -N <no_plink>             Include flag to skip generating plink separated files (default: set to generate split plink files). 
  -h                        show this help and exit

Example:
  bash prs_pipeline/src/run_split_plink_data.sh -1 /projects/standard/gdc/public/prs_methods/data/test/sim_1/AFR_simulation_study_sample -2 /projects/standard/gdc/public/prs_methods/data/test/sim_1/EUR_simulation_study_sample -t 40 -T 30 -v 30 -S 40

  Using default settings but skipping generation of separated files
    bash prs_pipeline/src/run_split_plink_data.sh -N
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
mkdir -p "${base_location}/samples_training"
mkdir -p "${base_location}/samples_testing"

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

  mv "${base_location}"/*samples.txt "${base_location}/samples_training"
  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed} \
    --no_plink
  mv "${base_location}"/*samples.txt "${base_location}/samples_testing"
else

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}
  mv "${base_location}"/*samples.txt "${base_location}/samples_training"

  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}
  mv "${base_location}"/*samples.txt "${base_location}/samples_testing"
fi

conda deactivate
conda deactivate


echo "Finished generating inputs with the following parameters  
training model: ${train_percent}%,
validation model: ${valid_percent}%, 
testing model: ${test_percent}% "

