input_path = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata.h5ad")

output_path = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class_ldata_dynamic.h5ad")

df <- data.frame(INPUT = input_path, OUTPUT = output_path)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/recover_dynamic/meta.csv",
    row.names = F)
