import scanpy as sc
import pandas as pd
full = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw.h5ad")
filtered = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered.h5ad")
adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")
df_full = full.obs.sampleid.value_counts().rename_axis('smapleID').reset_index(name='Number of Cells Sequenced')
df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table2/Number_of_Cells_Sample_part1.csv")
df.index = df["Unnamed: 0"].values
df = df.drop(["Unnamed: 0","samples"], axis=1)
df["smapleID"] = df.index
df = pd.merge(df_full, df, on='smapleID', how='outer')
temp = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/ATAC_filtered_cells/ATAC_Filtered_cell.csv",sep = " ")
temp.x = [x.replace("#", "_") for x in temp.x]
filtered = filtered[filtered.obs.index.isin(temp.x),]
df3 = filtered.obs.sampleid.value_counts().rename_axis('smapleID').reset_index(name='Pass Filter 3')
df = pd.merge(df,df3, on='smapleID', how='outer')
df4 = adata.obs.sampleid.value_counts().rename_axis('smapleID').reset_index(name='Pass Filter 4')
df = pd.merge(df,df4, on='smapleID', how='outer')

df.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table2/Number_of_Cells_Sample_part2.csv",index = False)