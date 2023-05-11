import loompy

loom_files = [
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_FR_2.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_11w2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_13W_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_13W_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_14w5d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_14w5d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_19W4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_19W4d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_20W2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_20W2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_23w1d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multi_Fetal_23w1d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_10w_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_10w_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_12w3d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_12w3d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_14w2d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_14w2d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_16w4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_16w4d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_20w1d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_20w1d_NR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_23w4d_FR.loom", 
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/Multiome_23w4d_NR.loom",
"/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/10x3_Lobe_19_D019_NeuN-v7.loom"
]
loompy.combine(
    loom_files,
    "/storage/singlecell/zz4/fetal_snakemake/results/RNA_velocity/combined.loom",
)