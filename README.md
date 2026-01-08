# prs_pipeline
A repository to aid in the use of multiple PRS models. 

Assumes the inputs are plink formatted and split by ancestry. Also assumes they are providing pc calculations as an input

Sequence to run scripts if using simulated data and needing to generate summary stats files
1. split_top_n_subjs.sh
2. run_split_plink_data.sh
3. generate_summary_stat_files.sh
4. restructure_output_dir.sh

run_prepare_prs.sh script runs the four modules to prepare for running a prs method in sequence for you.

```{bash}
Sample call with a smaller data set
sbatch prs_pipeline/run_prepare_prs.sh -1 /scratch.global/baron063/simulations/temp/AFR_simulation -2 /scratch.global/baron063/simulations/temp/EUR_simulation -P /projects/standard/gdc/public/prs_methods/data/adjusted_1kgPCs.tsv -n 300
```


Work on diverse ancestry PRS pipeline (if they need to be combined do so in the prs method script)
* CTSLEB - plink files are together for anc1 and anc2, They are separated for summary statitics and reference panels. 
* PRScsx - plink files could be together or separate for anc1 and anc2. Summary statistic files should be separate.
* prosper - plink files are together for anc1 and anc2
* TL-PRS - plink files are separate 
* viprs
