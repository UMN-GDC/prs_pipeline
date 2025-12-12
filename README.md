# prs_pipeline
A repository to aid in the use of multiple PRS models. 

Assumes the inputs are plink formatted and split by ancestry. Also assumes they are providing pc calculations as an input

Sequence to run scripts if using simulated data and needing to generate summary stats files
1. split_top_n_subjs.sh
2. run_split_plink_data.sh
3. generate_summary_stat_files.sh
4. restructure_output_dir.sh

