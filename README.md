<details open>
  <summary>Data Meta</summary>
|Sample.ID             |Donor.ID|Time |Region    |Days|Data.Type |
|----------------------|--------|-----|----------|----|----------|
|sn_multiome_d59       |1       |8w3d |Whole Eye |59  |Multiomics|
|Multiome_10w_FR       |2       |10w  |Macula    |70  |Multiomics|
|Multiome_10w_NR       |2       |10w  |Peripheral|70  |Multiomics|
|sn_multiome_d76c      |3       |10w6d|Macula    |76  |Multiomics|
|sn_multiome_d76p      |3       |10w6d|Peripheral|76  |Multiomics|
|Multi_Fetal_11w2d_FR  |4       |11w2d|Macula    |79  |Multiomics|
|Multi_Fetal_11w2d_FR_2|4       |11w2d|Macula    |79  |Multiomics|
|Multi_Fetal_11w2d_NR  |4       |11w2d|Peripheral|79  |Multiomics|
|Multiome_12w3d_FR     |5       |12w3d|Macula    |87  |Multiomics|
|Multiome_12w3d_NR     |5       |12w3d|Peripheral|87  |Multiomics|
|Multi_Fetal_13W_FR    |6       |13w  |Macula    |91  |Multiomics|
|Multi_Fetal_13W_NR    |6       |13w  |Peripheral|91  |Multiomics|
|Multiome_14w2d_FR     |7       |14w2d|Macula    |100 |Multiomics|
|Multiome_14w2d_NR     |7       |14w2d|Peripheral|100 |Multiomics|
|Multi_Fetal_14w5d_FR  |8       |14w5d|Macula    |103 |Multiomics|
|Multi_Fetal_14w5d_NR  |8       |14w5d|Peripheral|103 |Multiomics|
|Multiome_16w4d_FR     |9       |16w4d|Macula    |116 |Multiomics|
|Multiome_16w4d_NR     |9       |16w4d|Peripheral|116 |Multiomics|
|Multi_Fetal_19W4d_FR  |10      |19w4d|Macula    |137 |Multiomics|
|Multi_Fetal_19W4d_NR  |10      |19w4d|Peripheral|137 |Multiomics|
|Multiome_20w1d_FR     |11      |20w1d|Macula    |141 |Multiomics|
|Multiome_20w1d_NR     |11      |20w1d|Peripheral|141 |Multiomics|
|Multi_Fetal_20W2d_FR  |12      |20w2d|Macula    |142 |Multiomics|
|Multi_Fetal_20W2d_NR  |12      |20w2d|Peripheral|142 |Multiomics|
|Multi_Fetal_23w1d_FR  |13      |23w1d|Macula    |162 |Multiomics|
|Multi_Fetal_23w1d_NR  |13      |23w1d|Peripheral|162 |Multiomics|
|Multiome_23w4d_FR     |14      |23w4d|Macula    |165 |Multiomics|
|Multiome_23w4d_NR     |14      |23w4d|Peripheral|165 |Multiomics|
</details>

#################################################################################################### Annotation 
## Run Seurat QC 
## Apply QC cutoff to get a subset of cells 
## Run DoubletFinder on RNA-seq 
## Export DoubletFinder results to a csv file 
## Merge all samples and save as seurat object. 
## Merge all samples and save as anndata object, save the object and filtered object 
## Annote major cell types with adult data. Progenitors were labeled as MG. 
## Run scvi umap to annotate MG cells. Filter cells based on ATAC and RNA seq 
## Run subclass annotation within each major class. 
## Manually annotate subclass and use csv as input to update the object. 
## Oranize annotation by renaming the columns in obs. 
## Merge organized annotation 
## Run UMAP to check annotation 
## Save subclass annotation within each major class. 
## Run subclass annotation within each major class. 