import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Set the global font size
matplotlib.rcParams.update({'font.size': 20})  # Change 12 to your desired font size

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")
adata = adata[adata.obs.Region!="Whole Eye"]
adata = adata[adata.obs.majorclass=="Cone"]
sc.pp.normalize_total(adata,target_sum=1e6)
sc.pp.log1p(adata)

adata.obs["Weeks"] = adata.obs.Days.map(
    {
        59: "PCW8",
        70: "PCW10",
        76: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW15",
        103: "PCW15",
        116: "PCW15",
        137: "PCW19",
        141: "PCW19",
        142: "PCW19",
        162: "PCW23",
        165: "PCW23",
    }
)
markers = ['CNGB3', 'PDE6H']
df = sc.get.obs_df(adata, markers + ["Region","Days"])
grouped_df = df.groupby(["Region","Days"])
result = grouped_df.mean().reset_index()
result.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/plot_temporal_trend.csv")
df = result.loc[result.Region=="Macula",:]
# Convert the 'Date' column to datetime format

# Set the 'Date' column as the index
df.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.plot(df.index, df['CNGB3'], marker='o', linestyle='-', label='CNGB3')
plt.plot(df.index, df['PDE6H'], marker='s', linestyle='-', label='PDE6H')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 7)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/compare_pre_knowlege_Macula.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

df = result.loc[result.Region=="Peripheral",:]
# Convert the 'Date' column to datetime format

# Set the 'Date' column as the index
df.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.plot(df.index, df['CNGB3'], marker='o', linestyle='-', label='CNGB3')
plt.plot(df.index, df['PDE6H'], marker='s', linestyle='-', label='PDE6H')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 7)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure7/compare_pre_knowlege_Peripheral.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)