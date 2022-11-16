input_path = "/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object"
df <- data.frame("Samples_ID" = substr(list.files(input_path),1,nchar(list.files(input_path))-33),
                 "Samples" = list.files(input_path,full.names = T))

write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/QC_seurat_apply_filter/meta.csv",
          row.names = F)
