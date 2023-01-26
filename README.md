# Data 
|Samples               |Time |Region    |Days|Data Type |
|----------------------|-----|----------|----|----------|
|Multi_Fetal_11w2d_FR  |11w2d|Macula    |79  |Multiomics|
|Multi_Fetal_11w2d_FR_2|11w2d|Macula    |79  |Multiomics|
|Multi_Fetal_11w2d_NR  |11w2d|Peripheral|79  |Multiomics|
|Multi_Fetal_13W_FR    |13w  |Macula    |91  |Multiomics|
|Multi_Fetal_13W_NR    |13w  |Peripheral|91  |Multiomics|
|Multi_Fetal_14w5d_FR  |14w5d|Macula    |103 |Multiomics|
|Multi_Fetal_14w5d_NR  |14w5d|Peripheral|103 |Multiomics|
|Multi_Fetal_19W4d_FR  |19w4d|Macula    |137 |Multiomics|
|Multi_Fetal_19W4d_NR  |19w4d|Peripheral|137 |Multiomics|
|Multi_Fetal_20W2d_FR  |20w2d|Macula    |142 |Multiomics|
|Multi_Fetal_20W2d_NR  |20w2d|Peripheral|142 |Multiomics|
|Multi_Fetal_23w1d_FR  |23w1d|Macula    |162 |Multiomics|
|Multi_Fetal_23w1d_NR  |23w1d|Peripheral|162 |Multiomics|
|Multiome_10w_FR       |10w  |Macula    |70  |Multiomics|
|Multiome_10w_NR       |10w  |Peripheral|70  |Multiomics|
|Multiome_12w3d_FR     |12w3d|Macula    |87  |Multiomics|
|Multiome_12w3d_NR     |12w3d|Peripheral|87  |Multiomics|
|Multiome_14w2d_FR     |14w2d|Macula    |100 |Multiomics|
|Multiome_14w2d_NR     |14w2d|Peripheral|100 |Multiomics|
|Multiome_16w4d_FR     |16w4d|Macula    |116 |Multiomics|
|Multiome_16w4d_NR     |16w4d|Peripheral|116 |Multiomics|
|Multiome_20w1d_FR     |20w1d|Macula    |141 |Multiomics|
|Multiome_20w1d_NR     |20w1d|Peripheral|141 |Multiomics|
|Multiome_23w4d_FR     |23w4d|Macula    |165 |Multiomics|
|Multiome_23w4d_NR     |23w4d|Peripheral|165 |Multiomics|
|17W1D_Temporal_retina |17w1d|Temporal  |120 |snRNA-seq |
|17W1D_Nasal_retina    |17w1d|Nasal     |120 |snRNA-seq |
|17W1D_Fovea_retina    |17w1d|Macula    |120 |snRNA-seq |
|17w1d_I_Ret           |17w1d|Inferior  |120 |snRNA-seq |
|17w1d_S_Ret           |17w1d|Superior  |120 |snRNA-seq |
|multi_19W3D_I_ret     |19w3d|Inferior  |136 |Multiomics|
|multi_19W3D_N_RET     |19w3d|Nasal     |136 |Multiomics|
|multi_19W3d_T_ret     |19w3d|Temporal  |136 |Multiomics|
|multi_19w3d_F_ret     |19w3d|Macula    |136 |Multiomics|
|multi_19w3d_S_ret     |19w3d|Superior  |136 |Multiomics|

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
| QC_seurat_apply_filter | Apply filters based on QC results from QC_seurat.|
| DoubletFinder | Find Doublets. |
| DoubletFinder_filtered_cells | Save all filtered cell IDs into one csv file.|
| merge_h5ad | Merge all h5ad, add meta info, save, and save all filtered cells in another object.|
| cell_type_annotation_major_class |Run first round of major class annotation.| 
| run_umap_all_samples |Run second step of major class annotation.| 
#| scpred | Run scPredict() on filtered cells to infer cell type. |
#| merge_scpred_meta | Merge all meta data. |
#| run_umap_all_samples | Merge samples and run umap on all merged/filtered cells. |
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
|multivelo_recover_dynamics_run_umap| Run UMAP.|
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
| plot_umap_adult_annotated_object| Plot UMAP color by attributes for each cell type. |
| plot_umap_seperate_time| Plot UMAP color by time. |
| run_umap_seprate_by_cell_type_with_NRPC | Plot UMAP for each cell type branch. |



