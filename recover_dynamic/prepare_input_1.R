input_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap_ldata/annotated_umap.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata/merged_h5ad_adult_annotated_umap.h5ad"
)
output_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap_ldata_dynamic/annotated_umap_dynamic.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated_umap_ldata_dynamic/merged_h5ad_adult_annotated_umap.h5ad"
)

df <- data.frame("INPUT"=input_path,
                 "OUTPUT"=output_path)
write.csv(df,"/storage/singlecell/zz4/fetal_bash/scripts/recover_dynamic/meta.csv",row.names = F)
