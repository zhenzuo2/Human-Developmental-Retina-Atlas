input_path="/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"
samples <- list.dirs(input_path, full.names = F, recursive = F)
inputFiles <- paste(input_path,samples, '/outs/atac_fragments.tsv.gz', sep = "")
df <- data.frame("Samples"=samples,"inputFiles" = inputFiles)
df <- df[file.exists(df$inputFiles),]
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/ArchR_create_object/meta.csv",row.names = F)