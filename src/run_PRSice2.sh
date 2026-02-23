#!/bin/bash
summary_stats_file=$1
study_sample=$2
binary_flag=$3
phenotype_info_file=$4
output_path=$5
path_repo=$6

final_output=${output_path}/prs_pipeline/PRSice2/prs_method
mkdir -p $final_output

# run PRSice-2 on simulated data--note that this is the temporary data and is unlikely to produce correct results (use this more for formatting purposes)
Rscript "${path_repo}/src/PRSice.R" \
    --prsice "${path_repo}/src/PRSice_linux" \
    --base ${summary_stats_file} \
    --target ${study_sample} \
    --binary-target ${binary_flag} \
    --pheno ${phenotype_info_file} \
    --stat beta \
    --beta \
    --out ${final_output}/PRSice2_outputs
