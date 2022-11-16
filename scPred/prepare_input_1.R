input_path = "/storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/"
df <- data.frame("Samples_ID" = substr(list.files(input_path),1,nchar(list.files(input_path))-4),
                 "Samples" = list.files(input_path,full.names = T))
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/scPred/meta.csv",
          row.names = F)