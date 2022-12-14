sink("/storage/singlecell/zz4/fetal_bash/scripts/run_umap_seprate_by_cell_type_with_NRPC/run_umap_seprate_by_cell_type_with_NRPC.sh")
df <- read.csv("/storage/singlecell/zz4/fetal_bash/scripts/run_umap_seprate_by_cell_type_with_NRPC/meta.csv")
for (i in 1:nrow(df)) {
    cat(paste("slurmtaco.sh -p gpu -m 20G -t 1 -- python3 run_umap_seprate_by_cell_type_with_NRPC.py ",
        df$INPUT[i], " ", df$OUTPUT[i], " ", df$LABEL, ";\n", sep = ""))
}
sink()