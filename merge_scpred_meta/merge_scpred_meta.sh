input_meta_path="/storage/singlecell/zz4/fetal_bash/scripts/merge_scpred_meta/meta.csv"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/scPred_meta/"

slurmtaco.sh -p short -m 20G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript merge_scpred_meta.R "$input_meta_path" "$output_results_path";