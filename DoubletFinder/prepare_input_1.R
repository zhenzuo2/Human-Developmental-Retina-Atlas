input_path = "/storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object/"
df <- data.frame("Samples_ID" = substr(list.files(input_path),1,nchar(list.files(input_path))-39),
                 "Samples" = list.files(input_path,full.names = T))
                 
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/DoubletFinder/meta.csv",
          row.names = F)
