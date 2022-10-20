adata_file=(
/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/Cone_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/Rod_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult_umap/annotated_umap.h5ad
/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult_umap/annotated_umap.h5ad
)

output_fig_path=(
/storage/singlecell/zz4/fetal_bash/figures/AC_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/BC_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/Cone_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/RGC_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/Rod_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/MG_annotation_adult_umap/
/storage/singlecell/zz4/fetal_bash/figures/HC_annotation_adult_umap/
)

for i in ${!adata_file[@]}
do
    slurmtaco.sh -p short -m 2G -t 1 -- python3 plot_umap_adult_annotated_object.py "${adata_file[i]}" "${output_fig_path[i]}";
done