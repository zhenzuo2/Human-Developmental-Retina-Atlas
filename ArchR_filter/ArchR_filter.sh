input_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"
input_meta="/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/merged_raw_filtered_umap_10000_major_sub_class.obs.csv"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/ArchR/"

slurmtaco.sh -p short -m 100G -t 1 -- Rscript ArchR_filter.R "$input_path" "$input_meta" "$output_results_path";