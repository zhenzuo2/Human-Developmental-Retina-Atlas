output_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac.rds"

slurmtaco.sh --g00 -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript merge_atac.R "$output_file";