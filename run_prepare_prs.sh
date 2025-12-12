#!/bin/bash 
#SBATCH --time=13:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5g
#SBATCH --tmp=10g
#SBATCH -p agsmall
#SBATCH -o run_prepare_prs.out
#SBATCH --job-name prep_prs_pipeline


echo "This pipeline assumes that the path where you would like the data to be prepared is in the same directory as the first plink_file provided."
echo "This pipeline ."

### --- CONFIGURABLE DEFAULTS -------------------------------------------------
### split_top_n
# Defaults (can be overridden via CLI)
# Mandatory
plink_file_anc1="/projects/standard/gdc/public/prs_methods/data/test/sim_2/AFR_simulation"                     # target ancestry
plink_file_anc2="/projects/standard/gdc/public/prs_methods/data/test/sim_2/EUR_simulation"                     # training ancestry
p_pca=/projects/standard/gdc/public/prs_methods/data/simulated_1000G/adjusted_1kgPCs.tsv
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located

# For customization
gwas_number=120000
train_percent=50
valid_percent=20
test_percent=30
rand_seed=42
no_plink=0 # Values greater than 0 skip generation of plink data splits

### derived variables
plink_file_anc1_study_sample="${plink_file_anc1}_study_sample"                     # target ancestry # Default name from split_top_n
plink_file_anc2_study_sample="${plink_file_anc2}_study_sample"                     # training ancestry # Default name from split_top_n
base_location=$(dirname ${plink_file_anc1})
anc1_basename=$(basename ${plink_file_anc1})
anc1_gwas_input="${anc1_basename}_gwas"
anc2_basename=$(basename ${plink_file_anc2})
anc2_gwas_input="${anc2_basename}_gwas"

### restructure
# see anc1_basename & anc2_basename

### --- Usage for full script to be built ------------------------------------
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <plink_file_anc1>      Full path to target ancestry plink files (default: ${plink_file_anc1})
  -2 <plink_file_anc2>      Full path to training ancestry plink files  (default: ${plink_file_anc2})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -g <gwas_number>         Percent of data for GWAS (default: 120000)
  -h                        show this help and exit

Example:
  bash prs_pipeline/src/split_top_n_subjs.sh -1 /projects/standard/gdc/public/prs_methods/data/test/sim_1/AFR_simulation -2 /projects/standard/gdc/public/prs_methods/data/test/sim_1/EUR_simulation -g 100000

  Using default settings but changing number of samples being considered for the gwas section
    bash prs_pipeline/src/split_top_n_subjs.sh -g 15000
  
  For a smaller dataset
    bash prs_pipeline/src/split_top_n_subjs.sh -1 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/AFR_simulation -2 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/EUR_simulation -g 300
EOF
  exit 1
}

### --- ArgParser for this script ----------------------
### --- PARSE ARGS ------------------------------------------------------------
while getopts ":1:2:r:g:h" opt; do
  case "$opt" in
    1) plink_file_anc1="$OPTARG" ;;
    2) plink_file_anc2="$OPTARG" ;;
    r) path_to_repo="$OPTARG" ;;
    g) gwas_number="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}

### --- Actual script where each script gets sbatch --wait 

### ========== run_split_plink_data.sh ==========
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <plink_file_anc1_study_sample>      Full path to target ancestry plink files (default: ${plink_file_anc1_study_sample})
  -2 <plink_file_anc2_study_sample>      Full path to training ancestry plink files  (default: ${plink_file_anc2_study_sample})
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
  
  Small dataset example
    bash run_split_plink_data.sh -1 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/AFR_simulation_study_sample -2 /projects/standard/gdc/public/prs_methods/data/simulated_1000G/EUR_simulation_study_sample
EOF
  exit 1
}


### ========== generate_summary_stat_files.sh ==========
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <anc1_gwas_input>      Target ancestry plink files gwas split (default: ${anc1_gwas_input})
  -2 <anc2_gwas_input>      Training ancestry plink files gwas split  (default: ${anc2_gwas_input})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -b <base_location>        Full path to workspace (default: ${base_location})
  -S <seed>                 Randomization seed (default:42)
  -p <pca>                 Full path pca file (default: ${p_pca})
  -h                        show this help and exit

Example:
  bash prs_pipeline/src/generate_summary_stat_files.sh -1 AFR_simulation_gwas -2 EUR_simulation_gwas -b /projects/standard/gdc/public/prs_methods/data/test/sim_1

  Using default settings but changing base_location of data
    bash prs_pipeline/src/generate_summary_stat_files.sh -b /projects/standard/gdc/public/prs_methods/data/test/sim_2
EOF
  exit 1
}


### ========== restructure_output_dir.sh ==========
usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -1 <anc1_prefix>      Full path to target ancestry plink files (default: ${anc1_prefix})
  -2 <anc2_prefix>      Full path to training ancestry plink files  (default: ${anc2_prefix})
  -r <repo_path>            Path to prs_pipeline repo (default: ${path_to_repo})
  -b <base_location>         Percent of data for GWAS (default: ${base_location})
  -h                        show this help and exit

Example:
  bash prs_pipeline/src/generate_summary_stat_files.sh -1 AFR_simulation_gwas -2 EUR_simulation_gwas -b /projects/standard/gdc/public/prs_methods/data/test/sim_1

  Using default settings but changing base_location of data
    bash prs_pipeline/src/generate_summary_stat_files.sh -b /projects/standard/gdc/public/prs_methods/data/test/sim_2
EOF
  exit 1
}

