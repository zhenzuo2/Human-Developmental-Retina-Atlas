import loompy

loom_files = [
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_FR_2/Multi_Fetal_11w2d_FR_2.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_FR/Multi_Fetal_11w2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_NR/Multi_Fetal_11w2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_13W_FR/Multi_Fetal_13W_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_13W_NR/Multi_Fetal_13W_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_14w5d_FR/Multi_Fetal_14w5d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_14w5d_NR/Multi_Fetal_14w5d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_19W4d_FR/Multi_Fetal_19W4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_19W4d_NR/Multi_Fetal_19W4d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_20W2d_FR/Multi_Fetal_20W2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_20W2d_NR/Multi_Fetal_20W2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_23w1d_FR/Multi_Fetal_23w1d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_23w1d_NR/Multi_Fetal_23w1d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_10w_FR/Multiome_10w_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_10w_NR/Multiome_10w_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_12w3d_FR/Multiome_12w3d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_12w3d_NR/Multiome_12w3d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_14w2d_FR/Multiome_14w2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_14w2d_NR/Multiome_14w2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_16w4d_FR/Multiome_16w4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_16w4d_NR/Multiome_16w4d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_20w1d_FR/Multiome_20w1d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_20w1d_NR/Multiome_20w1d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_23w4d_FR/Multiome_23w4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_23w4d_NR/Multiome_23w4d_NR.loom"
]
loompy.combine(
    loom_files,
    "/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/combined.loom",
)