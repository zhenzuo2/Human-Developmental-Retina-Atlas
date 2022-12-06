input_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated.h5ad"
input_meta_file="/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/merged_h5ad_adult_annotated_obs.csv"
output_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/"

slurmtaco.sh --g01 -m 20G -t 1 -- python3 run_umap_saved_merged_adult_annotated_object.py "${input_file_path}" "${input_meta_file}" "${output_file_path}";
