input_path="/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"
output_path="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/"
meta_file="/storage/singlecell/zz4/fetal_bash/results/scPred_meta/meta.csv"

slurmtaco.sh --g01 -m 20G -t 1 -- python3.9 run_umap_all_samples.py "$input_path" "$output_path" "$meta_file";