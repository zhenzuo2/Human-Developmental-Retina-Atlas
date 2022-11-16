SAMPLES = ["Multi_Fetal_11w2d_FR", "Multi_Fetal_11w2d_FR_2",
    "Multi_Fetal_11w2d_NR", "Multi_Fetal_13W_FR", "Multi_Fetal_13W_NR",
    "Multi_Fetal_14w5d_FR", "Multi_Fetal_14w5d_NR", "Multi_Fetal_19W4d_FR",
    "Multi_Fetal_19W4d_NR", "Multi_Fetal_20W2d_FR", "Multi_Fetal_20W2d_NR",
    "Multi_Fetal_23w1d_FR", "Multi_Fetal_23w1d_NR", "Multiome_10w_FR",
    "Multiome_10w_NR", "Multiome_12w3d_FR", "Multiome_12w3d_NR", "Multiome_14w2d_FR",
    "Multiome_14w2d_NR", "Multiome_16w4d_FR", "Multiome_16w4d_NR", "Multiome_20w1d_FR",
    "Multiome_20w1d_NR", "Multiome_23w4d_FR", "Multiome_23w4d_NR"]
output_figures_path="figures/qc_seurat_object/"
output_results_path="results/qc_seurat_object/"
min_cells="10"
min_features="200"
i=list(range(1,26))
rule all:
    input:
        expand("{output_results_path}{sample}_min_cell_{min_cells}_min_features_{min_features}.rds", sample=SAMPLES,output_results_path=output_results_path,min_cells=min_cells,min_features=min_features),
        expand("{output_figures_path}{sample}_log_seurat_QC.svg", sample=SAMPLES,output_figures_path=output_figures_path),
        expand("{output_figures_path}{sample}_seurat_QC.svg", sample=SAMPLES,output_figures_path=output_figures_path)
        
rule QC_seurat_prepare_input_1:
    output:
        "scripts/QC_seurat/meta.csv"
    shell:
        "/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript scripts/QC_seurat/prepare_input_1.R"

rule QC_seurat_prepare_input_2:
    input:
        "scripts/QC_seurat/meta.csv"
    output:
        "scripts/QC_seurat/QC_seurat.sh"
    shell:
        "/storage/chen/home/zz4/anaconda3/envs/r/bin/Rscript scripts/QC_seurat/prepare_input_2.R"

rule QC_seurat:
    input:
        "scripts/QC_seurat/QC_seurat.sh"
    output:
        expand("{output_results_path}{sample}_min_cell_{min_cells}_min_features_{min_features}.rds", sample=SAMPLES,output_results_path=output_results_path,min_cells=min_cells,min_features=min_features),
        expand("{output_figures_path}{sample}_log_seurat_QC.svg", sample=SAMPLES,output_figures_path=output_figures_path),
        expand("{output_figures_path}{sample}_seurat_QC.svg", sample=SAMPLES,output_figures_path=output_figures_path)
    shell:
        """
        bash -c '
        cd scripts/QC_seurat/
        sh QC_seurat.sh'
        """
