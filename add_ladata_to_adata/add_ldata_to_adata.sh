adata_file=(
/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_with_label_umap/annotated_umap.h5ad,
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_with_label_umap/annotated_umap.h5ad
)

ldata_file="/storage/singlecell/zz4/fetal_bash/results/RNA_velocity/combined.loom"
output_file=(
/storage/singlecell/zz4/fetal_bash/results/merged_h5ad_adult_annotated/merged_h5ad_adult_annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap_ldata.h5ad
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap_ldata.h5ad
)

for i in ${!adata_file[@]}
do
    slurmtaco.sh -p gpu -m 20G -t 1 -- python3 add_ldata_to_adata.py "${adata_file[i]}" "${ldata_file}" "${output_file[i]}";
done;