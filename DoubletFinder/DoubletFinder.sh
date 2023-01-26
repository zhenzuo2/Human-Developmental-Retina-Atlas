input_path=/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/
output_results_path=/storage/singlecell/zz4/fetal_bash/results/DoubletFinder_seurat_object/
output_figures_path=/storage/singlecell/zz4/fetal_bash/figures/DoubletFinder_UMAP/
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R 17W1D_Fovea_retina /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//17W1D_Fovea_retina_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R 17w1d_I_Ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//17w1d_I_Ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R 17W1D_Nasal_retina /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//17W1D_Nasal_retina_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R 17w1d_S_Ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//17w1d_S_Ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R 17W1D_Temporal_retina /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//17W1D_Temporal_retina_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R multi_19w3d_F_ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//multi_19w3d_F_ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R multi_19W3D_I_ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//multi_19W3D_I_ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R multi_19W3D_N_RET /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//multi_19W3D_N_RET_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R multi_19w3d_S_ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//multi_19w3d_S_ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R multi_19W3d_T_ret /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//multi_19W3d_T_ret_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_11w2d_FR_2 /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_11w2d_FR_2_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_11w2d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_11w2d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_11w2d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_11w2d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_13W_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_13W_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_13W_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_13W_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_14w5d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_14w5d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_14w5d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_14w5d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_19W4d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_19W4d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_19W4d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_19W4d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_20W2d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_20W2d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_20W2d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_20W2d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_23w1d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_23w1d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multi_Fetal_23w1d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multi_Fetal_23w1d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_10w_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_10w_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_10w_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_10w_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_12w3d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_12w3d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_12w3d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_12w3d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_14w2d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_14w2d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_14w2d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_14w2d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_16w4d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_16w4d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_16w4d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_16w4d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_20w1d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_20w1d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_20w1d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_20w1d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_23w4d_FR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_23w4d_FR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
slurmtaco.sh -p short -m 20G -t 1 -- Rscript DoubletFinder.R Multiome_23w4d_NR /storage/singlecell/zz4/fetal_bash/results/after_qc_seurat_object//Multiome_23w4d_NR_nFeature_RNA_500_5000_MT_5_fitered.rds $output_results_path $output_figures_path;
