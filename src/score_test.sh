#!/bin/bash
#
# Score test data using trained PRS models and evaluate performance
# Usage: score_test.sh [options]

set -eu

usage() {
  cat <<EOF
Usage: $0 [options]

Required:
  --test-bfile <prefix>     PLINK prefix for test data
  --train-out-dir <dir>     Pipeline methods output directory (e.g., output_path/prs_pipeline)
  --sumstats <file>         GWAS summary stats file (needed for C+T)
  --path-repo <dir>         Path to prs_pipeline repo

Optional:
  --pheno-file <file>       Phenotype file (FID IID phenotype). Default: .fam column 6
  --test-pca-file <file>    PCA eigenvec file for test samples. Default: PLINK --pca 6
  --binary-flag <T/F>       Binary phenotype flag for PRSice2 (default: F)
  --ran-ct <true/false>     Whether C+T was run (default: false)
  --ran-ldpred2 <true/false> Whether LDpred2 was run (default: false)
  --ran-lassosum2 <true/false> Whether lassosum2 was run (default: false)
  --ran-prsice2 <true/false> Whether PRSice2 was run (default: false)
  -h, --help                Show this help
EOF
  exit 1
}

# --- Parse args ---
test_bfile=""
train_out_dir=""
pheno_file=""
test_pca_file=""
sumstats=""
path_repo=""
binary_flag=F
ran_ct=false
ran_ldpred2=false
ran_lassosum2=false
ran_prsice2=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --test-bfile)     test_bfile="$2";     shift 2 ;;
    --train-out-dir)  train_out_dir="$2";  shift 2 ;;
    --pheno-file)     pheno_file="$2";     shift 2 ;;
    --test-pca-file)  test_pca_file="$2";  shift 2 ;;
    --sumstats)       sumstats="$2";       shift 2 ;;
    --path-repo)      path_repo="$2";      shift 2 ;;
    --binary-flag)    binary_flag="$2";    shift 2 ;;
    --ran-ct)         ran_ct="$2";         shift 2 ;;
    --ran-ldpred2)    ran_ldpred2="$2";    shift 2 ;;
    --ran-lassosum2)  ran_lassosum2="$2";  shift 2 ;;
    --ran-prsice2)    ran_prsice2="$2";    shift 2 ;;
    -h|--help)        usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$test_bfile" || -z "$train_out_dir" || -z "$sumstats" || -z "$path_repo" ]]; then
  echo "ERROR: --test-bfile, --train-out-dir, --sumstats, and --path-repo are required."
  usage
fi

eval_dir="${train_out_dir}/test_evaluation"
mkdir -p "$eval_dir"

echo "[score_test] Evaluating PRS models on test data: $test_bfile"
echo "[score_test] Output dir: $eval_dir"

# --- 1. PCA for test samples ---
if [[ -n "$test_pca_file" ]]; then
  echo "[score_test] Using user-provided PCA file: $test_pca_file"
  pca_file="$test_pca_file"
else
  echo "[score_test] Computing PCA on test sample..."
  plink --bfile "$test_bfile" --pca 6 --allow-no-sex --out "${eval_dir}/pca"
  pca_file="${eval_dir}/pca.eigenvec"
fi

# --- 2. Phenotype for test samples ---
if [[ -n "$pheno_file" ]]; then
  echo "[score_test] Using user-provided phenotype file: $pheno_file"
  pheno_input="$pheno_file"
else
  echo "[score_test] Extracting phenotype from .fam column 6"
  awk 'BEGIN{print "FID\tIID\tphenotype"} {print $1, $2, $6}' OFS="\t" "${test_bfile}.fam" > "${eval_dir}/pheno.txt"
  pheno_input="${eval_dir}/pheno.txt"
fi

# --- 3. Score + evaluate each trained method ---
evaluate() {
  local label="$1"
  local out_prefix="$2"
  local profile_file="${eval_dir}/${label}.profile"

  if [[ ! -f "$profile_file" ]]; then
    echo "[score_test] WARNING: Profile not found: $profile_file — skipping ${label}"
    return
  fi

  Rscript "${path_repo}/src/evaluate_test.R" \
    "$profile_file" \
    "$pheno_input" \
    "$pca_file" \
    "${eval_dir}/${label}_results.txt" \
    "${eval_dir}/${label}_scores.txt"
}

# --- C+T ---
if [[ "$ran_ct" == true ]]; then
  ct_dir="${train_out_dir}/CT"
  ct_results="${ct_dir}/CT_prs_results.txt"
  if [[ -f "$ct_results" ]]; then
    echo "[score_test] Evaluating C+T..."
    # Find best threshold (max R2)
    best_p=$(awk 'NR>1 && $2+0 > max {max=$2; best=$1} END {print best}' "$ct_results")
    if [[ -z "$best_p" ]]; then
      echo "[score_test] WARNING: Could not parse best p-value from $ct_results — skipping C+T"
    else
      echo "[score_test]   Best C+T p-value threshold: $best_p"
      # Create single-range file at best threshold
      echo "${best_p} 0 ${best_p}" > "${eval_dir}/CT_best_range.txt"
      # Recreate SNP.pvalue from sumstats
      awk 'NR==1 {print "SNP pvalue"; next} {print $1, $8}' OFS="\t" "$sumstats" > "${eval_dir}/CT_SNP.pvalue"
      # Extract clumped SNP list from training temp dir
      valid_snp="${ct_dir}/temp/temp.valid.snp"
      if [[ ! -f "$valid_snp" ]]; then
        echo "[score_test] WARNING: Clumped SNP list not found: $valid_snp — skipping C+T"
      else
        plink --bfile "$test_bfile" \
          --score "$sumstats" 1 4 6 header \
          --q-score-range "${eval_dir}/CT_best_range.txt" "${eval_dir}/CT_SNP.pvalue" \
          --extract "$valid_snp" \
          --allow-no-sex \
          --out "${eval_dir}/CT"
        # q-score-range names the profile CT.{range_name}.profile; normalize to CT.profile
        actual_ct_profile=$(ls "${eval_dir}"/CT.*.profile 2>/dev/null | head -1)
        if [[ -f "$actual_ct_profile" ]]; then
          cp "$actual_ct_profile" "${eval_dir}/CT.profile"
        fi
        evaluate "CT" "${eval_dir}/CT"
      fi
    fi
  else
    echo "[score_test] WARNING: C+T results not found: $ct_results — skipping"
  fi
fi

# --- LDpred2 ---
if [[ "$ran_ldpred2" == true ]]; then
  ldpred2_base="${train_out_dir}/LDpred2/prs_method"
  echo "[score_test] Evaluating LDpred2..."
  for model in inf grid; do
    weights="${ldpred2_base}_${model}_weights.txt"
    if [[ -f "$weights" ]]; then
      plink --bfile "$test_bfile" \
        --score "$weights" 1 2 3 header \
        --allow-no-sex \
        --out "${eval_dir}/LDpred2_${model}"
      evaluate "LDpred2_${model}" "${eval_dir}/LDpred2_${model}"
    else
      echo "[score_test] WARNING: LDpred2 ${model} weights not found: $weights — skipping"
    fi
  done
fi

# --- lassosum2 ---
if [[ "$ran_lassosum2" == true ]]; then
  lassosum2_base="${train_out_dir}/lassosum2/prs_method"
  echo "[score_test] Evaluating lassosum2..."
  weights="${lassosum2_base}_weights.txt"
  if [[ -f "$weights" ]]; then
    plink --bfile "$test_bfile" \
      --score "$weights" 1 2 3 header \
      --allow-no-sex \
      --out "${eval_dir}/lassosum2"
    evaluate "lassosum2" "${eval_dir}/lassosum2"
  else
    echo "[score_test] WARNING: lassosum2 weights not found: $weights — skipping"
  fi
fi

# --- PRSice2 ---
if [[ "$ran_prsice2" == true ]]; then
  echo "[score_test] Evaluating PRSice2 on test data..."
  prsice2_out="${eval_dir}/PRSice2"
  mkdir -p "$prsice2_out"
  PRSice \
    --base "$sumstats" \
    --target "$test_bfile" \
    --binary-target "$binary_flag" \
    --pheno "$pheno_input" \
    --beta \
    --stat beta \
    --out "${prsice2_out}/PRSice2_outputs" 2>/dev/null || \
    echo "[score_test] WARNING: PRSice2 failed or not found in PATH"
  echo "[score_test]   PRSice2 output: ${prsice2_out}/PRSice2_outputs.*"
fi

echo "[score_test] Done. Results in: $eval_dir"
