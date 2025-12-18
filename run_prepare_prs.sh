#!/bin/bash 
#SBATCH --time=13:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH -o run_prepare_prs.out
#SBATCH --job-name prep_prs_pipeline

set -e 

echo "This pipeline assumes that the path where you would like the data to be prepared is in the same directory as the first plink_file provided."
echo "This pipeline ."

### --- CONFIGURABLE DEFAULTS -------------------------------------------------
### split_top_n
# Defaults (can be overridden via CLI)
# Mandatory
plink_file_anc1="/projects/standard/gdc/public/prs_methods/data/test/sim_2/AFR_simulation"                     # target ancestry
plink_file_anc2="/projects/standard/gdc/public/prs_methods/data/test/sim_2/EUR_simulation"                     # training ancestry
p_pca=/projects/standard/gdc/public/prs_methods/data/adjusted_1kgPCs.tsv
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located

# For customization
gwas_number=120000
train_percent=50
valid_percent=20
test_percent=30
rand_seed=42
no_plink=0 # Values greater than 0 skip generation of plink data splits

### restructure
# see anc1_basename & anc2_basename

### --- Usage for full script to be built ------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <plink_file_anc1>      Full path to target ancestry plink files (default: ${plink_file_anc1})
  -2 <plink_file_anc2>      Full path to training ancestry plink files (default: ${plink_file_anc2})
  -P <P_pca>                Full path to the pca or covariate file for use in gwas (default: ${p_pca}).
  -R <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -n <gwas_number>          Amount of data for GWAS (default: 120000)
  -t <train_percent>        Training percent as an integer to split study data into (default 50)
  -v <valid_percent>        Validation percent as an integer to split study data into (default 20)
  -T <test_percent>         Testing percent as an integer to split study data into (default 30)
  -s <rand_seed>            Randomization seed (default 42)
  -N <no_plink>             Include flag to skip generating plink separated files (default: set to generate split plink files). 
  -h                        show this help and exit

Example:
  Using default settings but changing number of samples being considered for the gwas section
    bash prs_pipeline/run_prepare_prs.sh -n 100000

  For a smaller dataset
    bash prs_pipeline/run_prepare_prs.sh -1 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/AFR_simulation -2 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/EUR_simulation -g 300 -P /projects/standard/gdc/public/prs_methods/data/simulated_1000G/adjusted_1kgPCs.tsv
EOF
  exit 1
}

### --- ArgParser for this script ----------------------
### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:P:R:n:t:v:T:s:Nh" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    P) p_pca="$OPTARG" ;;
    R) path_to_repo="$OPTARG" ;;
    n) gwas_number="$OPTARG" ;;
    t) train_percent="$OPTARG" ;;
    v) valid_percent="$OPTARG" ;;
    T) test_percent="$OPTARG" ;;
    s) rand_seed="$OPTARG" ;;
    N) no_plink=1 ;;
    h) usage ;;
    *) usage ;;
  esac
done

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}


### derived variables
plink_file_anc1_study_sample="${plink_file_anc1}_study_sample"                     # target ancestry # Default name from split_top_n
plink_file_anc2_study_sample="${plink_file_anc2}_study_sample"                     # training ancestry # Default name from split_top_n
base_location=$(dirname ${plink_file_anc1})
anc1_basename=$(basename ${plink_file_anc1})
anc1_gwas_input="${anc1_basename}_gwas"
anc2_basename=$(basename ${plink_file_anc2})
anc2_gwas_input="${anc2_basename}_gwas"

if [[ ! -f "${p_pca}" ]]; then
  echo "Error: The file ${p_pca} is incorrect"
fi

### --- Actual script where each script gets sbatch --wait 
log "Starting split_top_n_subjs.sh"
sbatch --wait "${path_to_repo}"/src/split_top_n_subjs.sh -1 "${plink_file_anc1}" -2 "${plink_file_anc2}" -r "${path_to_repo}" -g "${gwas_number}"

log "Running run_split_plink_data.sh"
if [ ${no_plink} -gt 0 ]; then
  sbatch --wait "${path_to_repo}"/src/run_split_plink_data.sh -1 "${plink_file_anc1_study_sample}" -2 "${plink_file_anc2_study_sample}" -r "${path_to_repo}" -t "${train_percent}" -v "${valid_percent}" -T "${test_percent}" -S "${rand_seed}" -N
else
  sbatch --wait "${path_to_repo}"/src/run_split_plink_data.sh -1 "${plink_file_anc1_study_sample}" -2 "${plink_file_anc2_study_sample}" -r "${path_to_repo}" -t "${train_percent}" -v "${valid_percent}" -T "${test_percent}" -S "${rand_seed}"
fi

log "Running generate_summary_stat_files.sh"
sbatch --wait "${path_to_repo}"/src/generate_summary_stat_files.sh -1 "${anc1_gwas_input}" -2 "${anc2_gwas_input}" -r "${path_to_repo}" -b "${base_location}" -S "${rand_seed}" -p "${p_pca}"

log "Running restructure_output_dir.sh"
bash "${path_to_repo}/src/restructure_output_dir.sh" -1 "${anc1_basename}" -2 "${anc2_basename}" -r "${path_to_repo}" -b "${base_location}"

log "Files are ready for a PRS method"
