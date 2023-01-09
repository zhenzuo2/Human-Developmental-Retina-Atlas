input_path = c("/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/AC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/all.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/BC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/Cone.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/HC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/MG.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/NRPC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/PRPC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/RGC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/Rod.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/RPC.h5ad")

output_path = c("/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/AC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/all.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/BC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/Cone.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/HC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/MG.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/NRPC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/PRPC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/RGC.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/Rod.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata_dynamic/RPC.h5ad")

df <- data.frame(INPUT = input_path, OUTPUT = output_path)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/recover_dynamic/meta.csv",
    row.names = F)
