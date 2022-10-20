output_file_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/"

slurmtaco.sh -p short -m 20G -t 1 -- python3 merge_annotated_h5ad.py "$output_file_path";