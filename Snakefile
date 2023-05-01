import pandas as pd
meta=pd.read_csv("/storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/QC_seurat/meta.csv")
samples=meta['Samples_ID'].tolist()
sample2file=dict(zip(meta['Samples_ID'].tolist(), meta['Samples'].tolist()))

wildcard_constraints:
    sample="[^/]+"

def get_rawh5(wildcards):
    return sample2file[wildcards.sample] 

rule all:
    ## Put input file here to run corresponding rules
    input:
        "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw.h5ad",
        "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered.h5ad"
            

####################################################################################################
# Annotation
rule QC_seurat:
    ## Seurat QC
    input:
        get_rawh5
    output:
        result="/storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object/{sample}.rds",
        figure="/storage/singlecell/zz4/fetal_snakemake/figures/after_qc_seurat_object/{sample}_seurat_QC.svg",
    shell:
        """
        bash -c '
        source /storage/singlecell/zz4/mambaforge/etc/profile.d/conda.sh && conda activate r
        Rscript /storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/QC_seurat/QC_seurat.R {wildcards.sample} {input} 200 10 /storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object/ /storage/singlecell/zz4/fetal_snakemake/figures/after_qc_seurat_object/'
        """

rule QC_seurat_apply_filter:
    ## Apply QC cutoff to get a subset of cells
    input:
        "/storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object/{sample}.rds"
    output:
        result="/storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object_apply_filter/{sample}.rds",
        figure="/storage/singlecell/zz4/fetal_snakemake/figures/after_qc_seurat_object_apply_filter/{sample}_after_filtered_seurat_QC.svg",
    shell:
        """
        bash -c '
        source /storage/singlecell/zz4/mambaforge/etc/profile.d/conda.sh && conda activate r
        Rscript /storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/QC_seurat_apply_filter/QC_seurat_apply_filter.R {wildcards.sample} {input} /storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object_apply_filter/ /storage/singlecell/zz4/fetal_snakemake/figures/after_qc_seurat_object_apply_filter/'
        """

rule DoubletFinder:
    ## Run DoubletFinder
    input:
        "/storage/singlecell/zz4/fetal_snakemake/results/after_qc_seurat_object_apply_filter/{sample}.rds"
    output:
        result="/storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_seurat_object/{sample}.rds",
        figure="/storage/singlecell/zz4/fetal_snakemake/figures/DoubletFinder_UMAP/{sample}.svg",
    shell:
        """
        bash -c '
        source /storage/singlecell/zz4/mambaforge/etc/profile.d/conda.sh && conda activate r
        Rscript /storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/DoubletFinder/DoubletFinder.R {wildcards.sample} {input} /storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_seurat_object/ /storage/singlecell/zz4/fetal_snakemake/figures/DoubletFinder_UMAP/'
        """

rule DoubletFinder_save_filtered_cells:
    ## Export DoubletFinder results to a csv file
    input:
        expand(
            [
            "/storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_seurat_object/{sample}.rds",
            ],
            sample=samples
            )
    output:
        "/storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_filtered_cells/DoubletFinder_filtered_cells.csv"
    shell:
        """
        bash -c '
        source /storage/singlecell/zz4/mambaforge/etc/profile.d/conda.sh && conda activate r
        Rscript /storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/DoubletFinder_save_filtered_cells/DoubletFinder_save_filtered_cells.R'
        """

rule merge_feature_bc_matrix_to_h5ad:
    ## Merge all samples and save as anndata object, save the object and filtered object
    input:
        "/storage/singlecell/zz4/fetal_snakemake/results/DoubletFinder_filtered_cells/DoubletFinder_filtered_cells.csv"
    output:
        "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw.h5ad",
        "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered.h5ad"
    shell:
        """
        bash -c '
        source /storage/singlecell/zz4/mambaforge/etc/profile.d/conda.sh && conda activate python
        python3.9 /storage/singlecell/zz4/fetal_snakemake/scripts/Annotation/merge_feature_bc_matrix_to_h5ad/merge_feature_bc_matrix_to_h5ad.py'
        """
