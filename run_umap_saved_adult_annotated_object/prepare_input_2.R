sink("/storage/singlecell/zz4/fetal_bash/scripts/run_umap_saved_adult_annotated_object/run_umap_saved_adult_annotated_object.sh")
df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/run_umap_saved_adult_annotated_object/meta.csv")
for (i in 1:nrow(df)) {
    cat(paste("slurmtaco.sh -p gpu -m 20G -t 1 -- python3 run_umap_saved_adult_annotated_object.py ",
        df$INPUT[i], " ", df$OUTPUT[i], ";\n", sep = ""))
}
sink()