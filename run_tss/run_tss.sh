input_atac_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac.rds"
output_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss.rds"
slurmtaco.sh --g01 -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript run_tss.R  ${input_atac_file} ${output_file};
