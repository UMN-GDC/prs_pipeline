#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30GB
#SBATCH --time=48:00:00
#SBATCH -p msismall
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=baron063@umn.edu 
#SBATCH -o VIprs_sims.out
#SBATCH --job-name VIprs_sims

source /projects/standard/gdc/public/envs/load_miniconda3.sh 
module load plink


####### CONFIG #######

path_data=/projects/standard/gdc/public/prs_methods/outputs/simulation_AFR_target_EUR_training_simPheno
targ_ancestry=AFR

# Inputs
piece1=${path_data}/target_sumstats_prosper_1_CTSLEB.txt      # CHR SNP BP A1 beta beta_se p rs_id
piece2=${path_data}/target_sumstats_prosper_1.txt             # rsid chr a1 a0 beta beta_se n_eff
piece3=${path_data}/target_sumstats_corrected_1.txt.assoc.linear # CHR SNP BP A1 TEST NMISS BETA STAT P

# Output
out_file=${path_data}/target_sumstats_viprs_1.txt

######################
# STEP 0: PLINK allele frequencies
######################
plink --bfile ${path_data}/${targ_ancestry}_sim_pheno.QC8 \
      --freq --out ${path_data}/${targ_ancestry}_sim_pheno_QC8_1_freq

# Extract SNP and frequency (trim spaces)
awk 'NR>1 {gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2,$5}' \
    ${path_data}/${targ_ancestry}_sim_pheno_QC8_1_freq.frq > /tmp/freqs.clean

######################
# STEP 1: Extract POS, T_STAT, P from piece3
######################
# piece3 columns: CHR SNP BP A1 TEST NMISS BETA STAT P
# Extract SNP, CHR, BP, STAT, P
awk 'NR>1 {print $2,$1,$3,$8,$9}' ${piece3} > /tmp/pos_tstat_pval.clean
# Columns: SNP CHR POS T_STAT P

######################
# STEP 2: Clean piece2 (trim whitespace)
######################
awk '{gsub(/^[ \t]+|[ \t]+$/,""); print}' ${piece2} > /tmp/piece2.clean

######################
# STEP 3: Merge all into VIPRS format
######################
awk 'BEGIN {OFS="\t"}
     NR==FNR {
         # Read POS, CHR, T_STAT, P
         snp=$1; chr_pos[$1]=$2; pos[$1]=$3; tstat[$1]=$4; pval[$1]=$5
         next
     }
     FILENAME==ARGV[2] {
         # Read A1_FREQ from PLINK
         snp=$1; freq[$1]=$2
         next
     }
     FNR==1 {
         # Print header
         print "#CHROM","POS","ID","REF","ALT1","A1","A1_FREQ","OBS_CT","BETA","SE","T_STAT","P"
         next
     }
     {
         snp=$1
         chr_val=$2
         a1=$3
         a0=$4
         beta=$5
         se=$6
         n=$7

         # Lookup POS, CHR from piece3
         pos_val = (snp in pos ? pos[snp] : "NA")
         chr_val2= (snp in chr_pos ? chr_pos[snp] : chr_val)

         # Lookup T_STAT, P from piece3
         t_val   = (snp in tstat ? tstat[snp] : "NA")
         p_val   = (snp in pval ? pval[snp] : "NA")

         # Lookup allele frequency
         freq_val= (snp in freq ? freq[snp] : "NA")

         print chr_val2, pos_val, snp, a0, a1, a1, freq_val, n, beta, se, t_val, p_val
     }' /tmp/pos_tstat_pval.clean /tmp/freqs.clean /tmp/piece2.clean > ${out_file}

######################
# STEP 4: Clean up temporary files
######################
rm -f /tmp/freqs.clean /tmp/piece2.clean /tmp/pos_tstat_pval.clean

echo "VIPRS-ready summary statistics saved to ${out_file}"


######################
# STEP 5: Splitting by chr to match example data
######################
# Path to your VIPRS-ready summary stats
in_file=${path_data}/target_sumstats_viprs_1.txt
out_dir=${path_data}/viprs_by_chr

mkdir -p $out_dir

# Get the header line
header=$(head -n1 $in_file)

# Use awk to split by chromosome
awk -v header="$header" '
NR>1 {          # skip header in input
    chr=$1
    print > "'"$out_dir"'/chr"chr".txt"
}
' $in_file

# Prepend header to each chromosome file
for f in $out_dir/chr*.txt; do
    sed -i "1i$header" "$f"
done

provided_sumstats="$out_dir/chr*"

#### Failing below because of some wierd way the sumstats file was built... trying again ####

viprs_fit -l "/projects/standard/gdc/public/prs_methods/ref/ref_viprs/EUR/chr_*" \
  -s "${in_file}" \
  --sumstats-format "plink" \
  --output-dir "/projects/standard/gdc/public/prs_methods/outputs/viprs"

for chr in {1..22}; do
    sed -i '1s/^#//' /projects/standard/gdc/public/prs_methods/outputs/simulation_AFR_target_EUR_training_simPheno/viprs_by_chr/chr${chr}.txt
done


python stable_archive/viprs_evaluate_pseudo_r2.py \
  --sumstats_files "${provided_sumstats}" \
  --sumstats_format plink \
  --inferred_betas "/projects/standard/gdc/public/prs_methods/outputs/viprs/VIPRS_EM.fit.gz"
