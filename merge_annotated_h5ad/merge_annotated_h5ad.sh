output_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/"
slurmtaco.sh -p short -m 20G -t 1 -- python3 merge_annotated_h5ad.py "$output_file_path";