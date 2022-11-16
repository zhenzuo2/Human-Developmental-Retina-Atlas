adata_file = "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_annotated_umap.h5ad"

cell_type = c("AC", "BC", "RGC", "HC", "Cone", "RPC")

meta_cluster_adata_file = c("/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult/merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult/merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult/merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult/merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/merged_object.h5ad",
                            "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult/merged_object.h5ad")

subclass_reference_file = c("/storage/singlecell/zz4/fetal_bash/data/manual_reference/AC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/BC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/RGC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/HC_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/Cone_subclass_annotation.csv",
                            "/storage/singlecell/zz4/fetal_bash/data/manual_reference/MG_subclass_annotation.csv")
majorclass_reference_file = "/storage/singlecell/zz4/fetal_bash/data/cell_subclass_majorclass_mapping.csv"

output_file_path = c("/storage/singlecell/zz4/fetal_bash/results/AC_annotation_adult/",
                     "/storage/singlecell/zz4/fetal_bash/results/BC_annotation_adult/",
                     "/storage/singlecell/zz4/fetal_bash/results/RGC_annotation_adult/",
                     "/storage/singlecell/zz4/fetal_bash/results/HC_annotation_adult/",
                     "/storage/singlecell/zz4/fetal_bash/results/Photoreceptor_annotation_adult/",
                     "/storage/singlecell/zz4/fetal_bash/results/MG_annotation_adult/")

output_fig_path = c("/storage/singlecell/zz4/fetal_bash/figures/AC_annotation_adult/",
                    "/storage/singlecell/zz4/fetal_bash/figures/BC_annotation_adult/",
                    "/storage/singlecell/zz4/fetal_bash/figures/RGC_annotation_adult/",
                    "/storage/singlecell/zz4/fetal_bash/figures/HC_annotation_adult/",
                    "/storage/singlecell/zz4/fetal_bash/figures/Photoreceptor_annotation_adult/",
                    "/storage/singlecell/zz4/fetal_bash/figures/MG_annotation_adult/")

n = length(cell_type)
sink("/storage/singlecell/zz4/fetal_bash/scripts/add_adult_annotated_reference_to_object/add_adult_annotated_reference_to_object.sh")

for (i in 1:n) {
  cat(paste("slurmtaco.sh -p short -m 20G -t 1 -- python3 add_adult_annotated_reference_to_object.py ",
            adata_file, " ", cell_type[i], " ", meta_cluster_adata_file[i],
            " ", subclass_reference_file[i], " ", majorclass_reference_file,
            " ", output_file_path[i], " ", output_fig_path[i], "\n", sep = ""))
}
sink()
