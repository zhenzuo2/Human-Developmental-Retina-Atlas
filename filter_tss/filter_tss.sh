input_atac_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss.rds"
output_file="/storage/singlecell/zz4/fetal_bash/results/merged_atac/atac_tss_filtered.rds"
slurmtaco.sh -p gpu -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript filter_tss.R  ${input_atac_file} ${output_file};
