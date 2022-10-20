input_path="/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"
output_file="/storage/singlecell/zz4/fetal_bash/results/merged_rna/merged_rna.rds"

slurmtaco.sh --g00 -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript merge_rna.R "$input_path" "$output_file";