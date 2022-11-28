adata_file = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap/merged_h5ad_adult_annotated_umap_X_scANVI.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap/annotated_umap.h5ad")

ldata_file = "/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.loom"

output_file = c("/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata/merged_h5ad_adult_annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
    "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad")

df <- data.frame(INPUT = adata_file, ldata_file = ldata_file, OUTPUT = output_file)
write.csv(df, "/storage/singlecell/zz4/fetal_bash/scripts/add_ldata_to_adata/meta.csv",
    row.names = F)