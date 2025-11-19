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
plink_file_anc1="/home/gdc/public/prs_methods/data/test/sim_1/AFR_simulation"                     # target ancestry
plink_file_anc2="/home/gdc/public/prs_methods/data/test/sim_1/EUR_simulation"                     # training ancestry

gwas_percent=40
train_percent=30
valid_percent=10
test_percent=20
rand_seed=42
no_plink=0
path_to_repo=/home/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located
# path_R_packages=/home/gdc/public/Ref/R
# export R_LIBS_USER="$path_R_packages"


### --- FUNCTIONS ------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <plink_file_anc1>      Full path to target ancestry plink files (default: ${plink_file_anc1})
  -2 <plink_file_anc2>      Full path to training ancestry plink files  (default: ${plink_file_anc2})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -g <gwas_percent>         Percent of data for GWAS (default: 40)
  -t <train_percent>        Percent of data for training (default: 30)
  -T <testing_percent>      Percent of data for testing (default: 20)
  -v <validation_percent>   Percent of data for validation (default: 10)
  -S <seed>                 Randomization seed (default:42)
  -N <no_plink>             Include flag to skip generating plink separated files (default: set to generate split plink files). 
  -h                        show this help and exit

Example:
  bash prs_pipeline/src/run_split_plink_data.sh -1 /home/gdc/public/prs_methods/data/test/sim_1/AFR_simulation -2 /home/gdc/public/prs_methods/data/test/sim_1/EUR_simulation -g 40 -t 30 -T 20 -v 10 -S 40

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
while getopts ":1:2:r:g:t:T:v:S:Nh" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    g) gwas_percent="$OPTARG" ;;
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


# Split phenotype / plink files
source /home/gdc/public/envs/load_miniconda3.sh
log "Splitting plink data for "$plink_file_anc1""
if [ "${no_plink}" -eq 1 ]; then
  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
    --gwas "${gwas_percent}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed} \
    --no_plink

  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --gwas "${gwas_percent}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed} \
    --no_plink
else

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
    --gwas "${gwas_percent}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}

  log "Splitting plink data for "$plink_file_anc2""

  python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
    --gwas "${gwas_percent}" \
    --train "${train_percent}" \
    --val "${valid_percent}" \
    --test "${test_percent}" \
    --seed ${rand_seed}
fi

conda deactivate
conda deactivate


echo "Finished generating inputs with the following parameters 
gwas: ${gwas_percent}%, 
training model: ${train_percent}%,
validation model: ${valid_percent}%, 
testing model: ${test_percent}% "

