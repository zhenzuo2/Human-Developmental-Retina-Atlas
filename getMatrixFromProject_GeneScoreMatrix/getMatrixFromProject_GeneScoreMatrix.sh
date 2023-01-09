input_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/Gene_score_matrix/"
slurmtaco.sh --g00 -m 100G -t 1 -- Rscript getMatrixFromProject_GeneScoreMatrix.R "$input_path" "$output_results_path";