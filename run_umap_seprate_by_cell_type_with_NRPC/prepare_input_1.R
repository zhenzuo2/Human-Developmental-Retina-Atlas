metas = c("/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/merged_h5ad_adult_annotated_obs.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/MG.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_AC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_BC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Cone.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_HC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_RGC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC_Rod.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/NRPC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/PRPC.csv",
    "/storage/singlecell/zz4/fetal_bash/results/merged_annotation_adult_with_label/RPC.csv")

labels = c("all", "MG", "AC", "BC", "Cone", "HC", "RGC", "Rod", "NRPC",
    "PRPC", "RPC")

output_file_path = "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/"

df <- data.frame(INPUT = metas, LABEL = labels, OUTPUT = output_file_path)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/run_umap_seprate_by_cell_type_with_NRPC/meta.csv",
    row.names = F)