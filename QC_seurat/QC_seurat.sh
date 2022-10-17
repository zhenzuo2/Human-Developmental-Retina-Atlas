input_path="/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/"
sample_id=(
Multi_Fetal_11w2d_FR
Multi_Fetal_11w2d_FR_2
Multi_Fetal_11w2d_NR
Multi_Fetal_13W_FR
Multi_Fetal_13W_NR
Multi_Fetal_14w5d_FR
Multi_Fetal_14w5d_NR
Multi_Fetal_19W4d_FR
Multi_Fetal_19W4d_NR
Multi_Fetal_20W2d_FR
Multi_Fetal_20W2d_NR
Multi_Fetal_23w1d_FR
Multi_Fetal_23w1d_NR
Multiome_10w_FR
Multiome_10w_NR
Multiome_12w3d_FR
Multiome_12w3d_NR
Multiome_14w2d_FR
Multiome_14w2d_NR
Multiome_16w4d_FR
Multiome_16w4d_NR
Multiome_20w1d_FR
Multiome_20w1d_NR
Multiome_23w4d_FR
Multiome_23w4d_NR
17W1D_Fovea_retina
17W1D_Nasal_retina
17W1D_Temporal_retina
)
min_cells=3
min_features=200
output_figures_path="/storage/singlecell/zz4/fetal_bash/figures/qc_seurat_object/"
output_results_path="/storage/singlecell/zz4/fetal_bash/results/qc_seurat_object/"


for f in ${sample_id[@]}
do
    slurmtaco.sh -p short -m 10G -t 1 -- /storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript QC_seurat.R "$input_path" "$f" "$min_cells" "$min_features" "$output_figures_path" "$output_results_path";
done
