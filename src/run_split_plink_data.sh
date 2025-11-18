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
# Allowed ancestry labels
readonly VALID_ANCS=("AFR" "AMR" "EAS" "EUR" "SAS")

# Defaults (can be overridden via CLI)
anc1="AFR"                     # target ancestry
anc2="EUR"                     # training ancestry
nsims=10                       # number of simulations to process

gwas_percent=0.4
TRAIN_MODEL_PERCENT=0.3
VALID_MODEL_PERCENT=0.1
TEST_MODEL_PERCENT=0.2
rand_seed=42
path_to_repo=/home/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located
path_R_packages=/home/gdc/public/Ref/R
export R_LIBS_USER="$path_R_packages"


### --- FUNCTIONS ------------------------------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <anc1>            target ancestry (default: ${anc1})
  -2 <anc2>            training ancestry (default: ${anc2})
  -n <nsims>           number of simulations to process (default: ${nsims})
  -s <0|1>             skip simulation generation: 1=yes, 0=no (default: ${skip_generate_sims})
  -b <base_location>   base output directory (default auto: simulation_\${anc1}_target_\${anc2}_training under public path)
  -r <repo_path>       path to TLPRS_tools repo (default: ${path_to_repo})
  -o <tmp_output>      temporary output location for QC steps (default: ${path_to_output_to})
  -R <RHO>             SNP effect covariance (0-1) (default: 0.8)
  -m <MAF>             maf (0-1) (default: 0.05)
  -H <HERIT>           herit (0-1) (default: 0.4)
  -S <seed>            Randomization seed (default:42)
  -h                   show this help and exit

Example:
  $0 -1 AFR -2 EUR -n 5 -b /path/to/output -s 0
EOF
  exit 1
}

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}

assert_file_exists() {
  local f="$1"
  if [[ ! -e "$f" ]]; then
    echo "Required file not found: $f" >&2
    exit 3
  fi
}

### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:r:g:t:T:v:S:h" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    g) gwas_percent="$OPTARG" ;;
    t) TRAIN_MODEL_PERCENT="$OPTARG" ;;
    T) TEST_MODEL_PERCENT="$OPTARG" ;;
    v) VALID_MODEL_PERCENT="$OPTARG" ;;
    S) rand_seed="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done


# Load environments / modules
# module load R/4.4.0-openblas-rocky8
module load plink
module load plink/2.00-alpha-091019


# Split phenotype / plink files
log "Splitting target/training data for run ${counter}"

source /home/gdc/public/envs/load_miniconda3.sh

echo "Splitting ancestry 1 data"
python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc1}" \
  --gwas "${gwas_percent}" \
  --train "${TRAIN_MODEL_PERCENT}" \
  --val "${VALID_MODEL_PERCENT}" \
  --test "${TEST_MODEL_PERCENT}" \
  --seed ${rand_seed}

echo "Splitting ancestry 2 data"
python "${path_to_repo}/src/split_plink_samples.py" "${plink_file_anc2}" \
  --gwas "${gwas_percent}" \
  --train "${TRAIN_MODEL_PERCENT}" \
  --val "${VALID_MODEL_PERCENT}" \
  --test "${TEST_MODEL_PERCENT}" \
  --seed ${rand_seed}

conda deactivate
conda deactivate


echo "Finished generating inputs with the following parameters 
gwas %: ${gwas_percent}, 
training model %: ${TRAIN_MODEL_PERCENT},
validation model %: ${VALID_MODEL_PERCENT}, 

