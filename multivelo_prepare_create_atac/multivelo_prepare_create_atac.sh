inputs=(
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_11w2d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_11w2d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_13W_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_10w_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_11w2d_FR_2
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_10w_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_13W_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_14w5d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_14w5d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_19W4d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_19W4d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_20W2d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_20W2d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_23w1d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multi_Fetal_23w1d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_12w3d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_12w3d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_14w2d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_14w2d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_16w4d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_16w4d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_20w1d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_20w1d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_23w4d_FR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/Multiome_23w4d_NR
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/multi_19w3d_F_ret
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/multi_19W3D_I_ret
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/multi_19W3D_N_RET
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/multi_19w3d_S_ret
/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/multi_19W3d_T_ret
)

samples=(
Multi_Fetal_11w2d_FR
Multi_Fetal_11w2d_NR
Multi_Fetal_13W_FR
Multiome_10w_FR
Multi_Fetal_11w2d_FR_2
Multiome_10w_NR
Multi_Fetal_13W_NR
Multi_Fetal_14w5d_FR
Multi_Fetal_14w5d_NR
Multi_Fetal_19W4d_FR
Multi_Fetal_19W4d_NR
Multi_Fetal_20W2d_FR
Multi_Fetal_20W2d_NR
Multi_Fetal_23w1d_FR
Multi_Fetal_23w1d_NR
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
multi_19w3d_F_ret
multi_19W3D_I_ret
multi_19W3D_N_RET
multi_19w3d_S_ret
multi_19W3d_T_ret
)

output_result_path=/storage/singlecell/zz4/fetal_bash/results/ATAC_aggregate_peaks_10x/

output_figure_path=/storage/singlecell/zz4/fetal_bash/figures/ATAC_aggregate_peaks_10x/

for i in "${!inputs[@]}"; do
    slurmtaco.sh --g01 -m 20G -t 1 -- python3 multivelo_prepare_create_atac.py ${inputs[i]} ${samples[i]} ${output_result_path} ${output_figure_path};
done

