mkdir /storage/singlecell/zz4/fetal_bash/figures/Pando_grn/

seurat_pando_object=(
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_All_feature_selection_FALSE.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_All_feature_selection_TRUE.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Macula_feature_selection_FALSE.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Macula_feature_selection_TRUE.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Peripheral_feature_selection_FALSE.rds
/storage/singlecell/zz4/fetal_bash/results/BC_Pando/BC_seurat_object_Peripheral_feature_selection_TRUE.rds
)
time_cds=(
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/monocle3.rds
)
labels=(
BC_All_F
BC_All_T
BC_Macula_F
BC_Macula_T
BC_Peripheral_F
BC_Peripheral_T
)
output_dir=(
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
/storage/singlecell/zz4/fetal_bash/figures/Pando_grn/BC/
)

for i in "${!seurat_pando_object[@]}"
do
	slurmtaco.sh --g01 -m 50G -t 1 -- Rscript plot_grn.R  ${seurat_pando_object[i]} ${time_cds[i]} ${labels[i]} ${output_dir[i]};
done


