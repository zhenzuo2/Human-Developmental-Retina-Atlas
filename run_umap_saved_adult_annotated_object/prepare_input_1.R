input_file_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label/AC_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label/BC_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Cone_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label/Rod_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label/HC_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label/RGC_major_sub_class.h5ad",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label/RPC_major_sub_class.h5ad"
)

output_file_path=c(
  "/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap/",
  "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap/"
)

df <- data.frame("INPUT" = input_file_path, "OUTPUT" = output_file_path)
write.csv(df,"/storage/singlecell/zz4/fetal_bash/scripts/run_umap_saved_adult_annotated_object/meta.csv", row.names = F)