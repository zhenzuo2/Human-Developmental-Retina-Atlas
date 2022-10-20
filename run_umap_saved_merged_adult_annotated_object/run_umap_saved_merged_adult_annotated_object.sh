input_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated.h5ad"

output_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/"

slurmtaco.sh -p gpu -m 20G -t 1 -- python3 run_umap_saved_merged_adult_annotated_object.py "${input_file_path}" "${output_file_path}";
