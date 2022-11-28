input_file_path=(
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_10w_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_10w_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_11w2d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_11w2d_FR_2.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_11w2d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_14w5d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_14w5d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_13W_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_13W_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_20W2d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_12w3d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_19W4d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_19W4d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_23w1d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_20W2d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_12w3d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_14w2d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multi_Fetal_23w1d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_16w4d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_16w4d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_20w1d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_23w4d_FR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_23w4d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_20w1d_NR.h5ad
/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/Multiome_14w2d_NR.h5ad
)

samples=(
Multiome_10w_NR
Multiome_10w_FR
Multi_Fetal_11w2d_FR
Multi_Fetal_11w2d_FR_2
Multi_Fetal_11w2d_NR
Multi_Fetal_14w5d_NR
Multi_Fetal_14w5d_FR
Multi_Fetal_13W_FR
Multi_Fetal_13W_NR
Multi_Fetal_20W2d_FR
Multiome_12w3d_FR
Multi_Fetal_19W4d_FR
Multi_Fetal_19W4d_NR
Multi_Fetal_23w1d_NR
Multi_Fetal_20W2d_NR
Multiome_12w3d_NR
Multiome_14w2d_FR
Multi_Fetal_23w1d_FR
Multiome_16w4d_FR
Multiome_16w4d_NR
Multiome_20w1d_FR
Multiome_23w4d_FR
Multiome_23w4d_NR
Multiome_20w1d_NR
Multiome_14w2d_NR
)

output_result_path=/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x_noramlized/

for i in "${!input_file_path[@]}"; do
    slurmtaco.sh --g01 -m 20G -t 1 -- python3 multivelo_normalize_atac_TF.py ${input_file_path[i]} ${samples[i]} ${output_result_path};
done

