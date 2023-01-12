# ATAC-seq analysis
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| ArchR_create_object  | Create ArchR object for all data with ATAC-seq. |
| ArchR_filter  | Filter Doublets. |
| getMatrixFromProject_GeneScoreMatrix | Extract gene score matrix from ArchR project. |
| save_gene_score_as_h5ad | Save gene score matrix as h5ad file. |  

# Data annotation
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| QC_seurat  | First-round of QC. Read10X_h5() and saveRDS(). Generate QC figures to find optimal QC parameters for each sample. |
| QC_seurat_apply_filter  | Apply filters based on QC results from QC_seurat.  |
| DoubletFinder | Find Doublets. |
| scpred | Run scPredict() on filtered cells to infer cell type. |
| merge_scpred_meta | Merge all meta data. |
| run_umap_all_samples | Run umap on all merged/filtered cells. |
| annotation_with_adult | Assign subclass labels with adult data. |
| add_adult_annotated_reference_to_object | Add adult annotated reference to the object. |
| merge_annotated_h5ad | Merge annotated object to one anndata object. |
| run_umap_adult_annotated_object | Run umap on objects with annotations from adults. |
| plot_umap_adult_annotated_object | Plot umap on objects with annotations from adults. |

# Merge samples
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| merge_atac | Merge atac and save as rds. |
| run_tss | Caculate tss score for merged atac rds. |
| filter_tss | Filter peaks based on tss. |
| merge_rna | Merge RNA and save as rds. |

# Volocity Analysis
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| velocyto | Run velocyto to get loom files for each sample. |
| merge_looms | Merge all looms from velocyto results. |
| convert_loom_to_h5ad | Convert merged loom file to anndata object. |
| add_ldata_to_adata | Add merged loom h5ad object to anndata object. |
| compute_vk | Kernel which computes a transition matrix based on RNA velocity. |
| infer_fate | Infer fate for all cells and save predicted terminal states. |
| infer_fate_NRPC | Infer fate for NRPCs. |
| infer_fate_NRPC_save_meta | Plot well known markers for each NRPC fate to validate prediction. |
| recover_dynamic | Recover dynamics for each cell type branch. |
| plot_dynamic | Plot scv.pl.velocity_embedding_stream(). |


# multivelo
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| multivelo_prepare_create_atac | Preprocessing the ATAC counts. |
| multivelo_normalize_atac_TF | Enhancer peaks were aggregated with 10x annotated promoter peaks and then normalized. |
| multivelo_shared_barcodes_RNA_ATAC | Find common cells. |
| multivelo_seurat_wnn | Run WNN on common cells. |
| multivelo_knn_smooth_chrom | Smooth ATAC-seq data based on WNN results. |
| multivelo_recover_dynamics_chrom | Recover dynamics with chrom. |

# Pando
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| Pando_merge_object | Data preparation for pando. |
| run_pando | Run Pando to infer GRNs. |
| plot_grn_region | Plot GRNs. |
| plot_grn_region_diff | Plot regional difference GRN. |

# DE analysis
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| monocle3_de_analysis | Run DE analysis. |
| plot_Region_DE_heatmap | Plot DE heatmap for regional DE genes. |

# Cell Cycle
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| cell_cycle_infer | Infer cell cycle for RPCs. |

# Others
| Folder Name  | Functions and aims |
| ------------- | ------------- |
| plot_cell_type_composition | Plot cell type composition. |
| summary_sample_meta | Summary sample meta information. | 





