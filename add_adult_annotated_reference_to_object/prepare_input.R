adata_file = "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000.h5ad"

cell_type = c("AC", "BC", "RGC", "HC", "Cone", "Rod","MG")

meta_cluster_adata_file = c("/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/AC_annotation_adult/AC_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/BC_annotation_adult/BC_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/RGC_annotation_adult/RGC_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/HC_annotation_adult/HC_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Cone_annotation_adult/Cone_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Rod_annotation_adult/Rod_merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/MG_annotation_adult/MG_merged_object.h5ad")

subclass_reference_file = c("/storage/singlecell/zz4/fetal_bash/data/manual_reference/AC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/BC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/RGC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/HC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/Cone_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/Rod_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/MG_subclass_annotation.csv")

majorclass_reference_file = "/storage/singlecell/zz4/fetal_bash/data/cell_subclass_majorclass_mapping.csv"

output_file_path = c("/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/AC_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/BC_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/RGC_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/HC_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Cone_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/Rod_annotation_adult/",
                            "/storage/singlecell/zz4/fetal_bash/results/subclass_annotation/MG_annotation_adult/")

output_fig_path = c("/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/AC_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/BC_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/RGC_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/HC_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/Photoreceptor_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/Photoreceptor_annotation_adult_with_label/",
                    "/storage/singlecell/zz4/fetal_bash/figures/subclass_annotation/MG_annotation_adult_with_label/")

n = length(cell_type)
sink("/storage/singlecell/zz4/fetal_bash/scripts/add_adult_annotated_reference_to_object/add_adult_annotated_reference_to_object.sh")

for (i in 1:n) {
  cat(paste("slurmtaco.sh -p gpu -m 20G -t 1 -- python3 add_adult_annotated_reference_to_object.py ",
            adata_file, " ", cell_type[i], " ", meta_cluster_adata_file[i],
            " ", subclass_reference_file[i], " ", majorclass_reference_file,
            " ", output_file_path[i], " ", output_fig_path[i], "\n", sep = ""))
}
sink()
