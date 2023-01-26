input_meta_file="/storage/singlecell/zz4/fetal_bash/scripts/ArchR_create_object/meta.csv"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"

slurmtaco.sh -p short -m 20G -t 10 -- Rscript creat_object.R "$input_meta_file" "$output_results_path";