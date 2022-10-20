# fetal_project_bash

## QC_seurat
Read10X_h5() and saveRDS(). Generate QC figures.

## QC_seurat_apply_filter
Apply filters based on QC results from QC_seurat.

## DoubletFinder
Find Doublets.

## scpred
Run scPredict().

## merge_scpred_meta
Merge all meta data.

## run_umap
Run umap on all merged/filtered cells.

## ArchR_create_object
Create ArchR object

## ArchR_filter
Filter Doublets.

## velocyto
Run velocyto

## merge_looms
Merge all looms from velocyto results.

## annotate_with_adult
Assign subclass labels with adult data

## save_adult_annotated_object
Save objects with annotations from adults

## run_umap_adult_annotated_object
Run umap on objects with annotations from adults

## plot_umap_adult_annotated_object
Plot umap on objects with annotations from adults

## merge_atac
Merge atac and save as rds

## merge_rna
Merge RNA and save as rds
