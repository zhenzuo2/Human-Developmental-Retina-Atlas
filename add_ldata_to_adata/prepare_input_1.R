adata_file = c("/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/AC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/all.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/BC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/Cone.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/HC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/MG.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/NRPC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/PRPC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/RGC.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/Rod.h5ad", 
               "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap/RPC.h5ad"
)

ldata_file = "/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.loom"

output_file = c("/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/AC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/all.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/BC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/Cone.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/HC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/MG.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/NRPC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/PRPC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/RGC.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/Rod.h5ad", 
                "/storage/singlecell/zz4/fetal_bash/results/annotation_adult_with_label_NRPC_umap_ldata/RPC.h5ad"
)

df <- data.frame(INPUT = adata_file, ldata_file = ldata_file, OUTPUT = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/add_ldata_to_adata/meta.csv",
          row.names = F)