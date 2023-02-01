adata_file = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad")

ldata_file = "/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.loom"

output_file = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata.h5ad")

df <- data.frame(INPUT = adata_file, ldata_file = ldata_file, OUTPUT = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/add_ldata_to_adata/meta.csv",
          row.names = F)