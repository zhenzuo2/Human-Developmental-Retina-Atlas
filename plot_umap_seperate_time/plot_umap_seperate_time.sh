input_adata="/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap.h5ad"
output_fig_path="/storage/singlecell/zz4/fetal_bash/figures/merged_seperate_time_point/"

slurmtaco.sh -p short -m 2G -t 1 -- python3 plot_umap_seperate_time.py ${input_adata} ${output_fig_path};
