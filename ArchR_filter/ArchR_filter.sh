input_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"
input_meta="/storage/singlecell/zz4/fetal_bash/results/scPred_meta/meta.csv"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"

slurmtaco.sh -p short -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript ArchR_filter.R "$input_path" "$input_meta" "$output_results_path";