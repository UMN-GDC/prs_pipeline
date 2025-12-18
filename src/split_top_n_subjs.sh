#!/bin/bash 
#SBATCH --time=1:00:00
#SBATCH --ntasks=2
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH -p agsmall
#SBATCH -o split_top_n_subjs.out
#SBATCH --job-name prs_split_top_subjs

### --- CONFIGURABLE DEFAULTS -------------------------------------------------
# Defaults (can be overridden via CLI)
plink_file_anc1="/projects/standard/gdc/public/prs_methods/data/test/sim_2/AFR_simulation"                     # target ancestry
plink_file_anc2="/projects/standard/gdc/public/prs_methods/data/test/sim_2/EUR_simulation"                     # training ancestry
#/projects/standard/gdc/public/prs_methods/data/simulated_1000G/
path_to_repo=/projects/standard/gdc/public/prs_methods/scripts/prs_pipeline # place the code is located

gwas_number=120000
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

log() {
  local msg="$1"
  echo "[$(date '+%F %T')] $msg"
}



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

### --- ENVIRONMENT ---------------------------------------------------------
module load plink

### --- DERIVED VARIABLES ---------------------------------------------------
plink_anc1_fam="${plink_file_anc1}".fam
plink_anc2_fam="${plink_file_anc2}".fam

output_anc1_list="${plink_file_anc1}"_gwas_list.txt
output_anc2_list="${plink_file_anc2}"_gwas_list.txt

plink_anc1_gwas="${plink_file_anc1}"_gwas
plink_anc2_gwas="${plink_file_anc2}"_gwas

plink_anc1_new="${plink_file_anc1}"_study_sample
plink_anc2_new="${plink_file_anc2}"_study_sample

extraction_function() {
    local gwas_num="$1"
    local input_fam="$2"
    local output_text_list="$3"

    head -n "$gwas_num" ${input_fam} | awk '{print $1, $2}' > ${output_text_list}
}

plink_split_function() {
    local input="$1"
    local text_list="$2"
    local output_file="$3"
    local output_complement="$4"

    plink --bfile ${input} \
          --keep ${text_list} \
          --make-bed \
          --out ${output_file}

    plink --bfile ${input} \
         --remove ${text_list} \
         --make-bed \
         --out ${output_complement}
}


### --- ACTUAL SCRIPT -------------------------------------------------------

log "Starting the process for provided anc1"
extraction_function "${gwas_number}" "${plink_anc1_fam}" "${output_anc1_list}"
plink_split_function "${plink_file_anc1}" "${output_anc1_list}" "${plink_anc1_gwas}" "${plink_anc1_new}"

log "Starting the process for provided anc2"
extraction_function "${gwas_number}" "${plink_anc2_fam}" "${output_anc2_list}"
plink_split_function "${plink_file_anc2}" "${output_anc2_list}" "${plink_anc2_gwas}" "${plink_anc2_new}"

log "End of split_top_n_subjs.sh script"
