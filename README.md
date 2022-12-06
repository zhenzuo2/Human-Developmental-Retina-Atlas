# fetal_project_bash

# ATAC-seq analysis

## ArchR_create_object
Create ArchR object for all data with ATAC-seq.

## ArchR_filter
Filter Doublets.

# Data annotation

## QC_seurat
First-round of QC. Read10X_h5() and saveRDS(). Generate QC figures to find optimal QC parameters for each sample.

## QC_seurat_apply_filter
Apply filters based on QC results from QC_seurat.

## DoubletFinder
Find Doublets.

## scpred
Run scPredict() on filtered cells to infer cell type.

## merge_scpred_meta
Merge all meta data.

## run_umap_all_samples
Run umap on all merged/filtered cells.

## annotate_with_adult
Assign subclass labels with adult data

## add_adult_annotated_reference_to_object
add_adult_annotated_reference_to_object

## merge_annotated_h5ad
Merge annotated object to one anndata object

## run_umap_adult_annotated_object
Run umap on objects with annotations from adults

## plot_umap_adult_annotated_object
Plot umap on objects with annotations from adults

# Merge data
## merge_atac
Merge atac and save as rds

## merge_rna
Merge RNA and save as rds

# Volocity Analysis

## velocyto
Run velocyto

## merge_looms
Merge all looms from velocyto results.

## compute_vk

## infer_fate


## multivelo_prepare_create_atac
Preprocessing the ATAC counts

## multivelo_normalize_atac_TF
Enhancer peaks were aggregated with 10x annotated promoter peaks and then normalized.

## multivelo_shared_barcodes_RNA_ATAC
Find common cells.

## multivelo_seurat_wnn
Run WNN on common cells.

## multivelo_knn_smooth_chrom
Smooth ATAC-seq data based on WNN results.

## multivelo_recover_dynamics_chrom




