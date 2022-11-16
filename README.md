# fetal_project_bash

## ArchR_create_object
Create ArchR object for all data with ATAC-seq.

## ArchR_filter
Filter Doublets.

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

## velocyto
Run velocyto

## merge_looms
Merge all looms from velocyto results.

## annotate_with_adult
Assign subclass labels with adult data

## add_adult_annotated_reference_to_object
add_adult_annotated_reference_to_object

## run_umap_adult_annotated_object
Run umap on objects with annotations from adults

## plot_umap_adult_annotated_object
Plot umap on objects with annotations from adults

## merge_atac
Merge atac and save as rds

## merge_rna
Merge RNA and save as rds
