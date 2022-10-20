adata_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap.h5ad"
ldata_file="/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.loom"
output_file="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap_ldata.h5ad"

slurmtaco.sh -p gpu -m 20G -t 1 -- python3 add_ldata_to_adata.py "${adata_file}" "${ldata_file}" "${output_file}";
