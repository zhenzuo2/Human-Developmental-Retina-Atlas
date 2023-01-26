output_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac.rds"

slurmtaco.sh --g00 -m 100G -t 1 -- Rscript merge_atac.R "$output_file";