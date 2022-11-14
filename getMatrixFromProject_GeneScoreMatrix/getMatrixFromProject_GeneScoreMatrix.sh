input_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"
output_results_path="c"
slurmtaco.sh --g00 -m 100G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript getMatrixFromProject_GeneScoreMatrix.R "$input_path" "$output_results_path";