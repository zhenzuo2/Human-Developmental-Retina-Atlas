"""Dev Cell. Author manuscript; available in PMC 2021 May 18.
Published in final edited form as:
Dev Cell. 2020 May 18; 53(4): 473–491.e9.
Published online 2020 May 7. doi: 10.1016/j.devcel.2020.04.009
PMCID: PMC8015270
NIHMSID: NIHMS1629633
PMID: 32386599
Single-Cell Analysis of Human Retina Identifies Evolutionarily Conserved and Species-Specific Mechanisms Controlling Development"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Set the global font size
matplotlib.rcParams.update({'font.size': 20})  # Change 12 to your desired font size


adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")

adata = adata[adata.obs.majorclass.isin(["PRPC"])]
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

markers = ['CYP26A1', 'DIO2', 'CDKN1A', 'ANXA2', 'FRZB',"CCN2"]

df = sc.get.obs_df(adata, markers + ["Region"])

df = df.set_index('Region').stack().reset_index()
df.columns = ['Region', 'gene', 'value']


df['NonZero'] = df['value'] != 0

# Group by 'Region' and 'Gene', then calculate the percentage of non-zero values
result = df.groupby(['Region', 'gene'])['NonZero'].mean().reset_index()

df = pd.DataFrame(result)

# Create subplots for each gene
genes = df['gene'].unique()
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16, 10), sharey=True)

# Define colors for each region
region_colors = {"Whole Eye": "#ffdd05", "Macula" : "#F8766D", "Peripheral": "#00BFC4"}

# Iterate over genes and create bar plots
for i, gene in enumerate(genes):
    ax = axes.flatten()[i]
    gene_df = df[df['gene'] == gene]
    for region in gene_df['Region'].unique():
        region_df = gene_df[gene_df['Region'] == region]
        ax.bar(region_df['Region'], region_df['NonZero']*100, color=region_colors.get(region), label=region)
    ax.set_title(gene)
    ax.set_ylabel('Cells (Percent) in PRC')
    ax.set_xlabel('')
    ax.legend()

# Adjust layout for better spacing
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)
################

adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad")
adata = adata[adata.obs.Region!="Whole Eye"]
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
markers = ['ALDH1A1', 'ALDH1A2','ALDH1A3','TBX5','VAX2']
df = sc.get.obs_df(adata, markers + ["Region","Days"])
grouped_df = df.groupby(["Region","Days"])
result = grouped_df.mean().reset_index()
result.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege.csv")
df = result.loc[result.Region=="Macula",:]
# Convert the 'Date' column to datetime format

# Set the 'Date' column as the index
df.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8,4)
plt.plot(df.index, df['ALDH1A1'], marker='o', linestyle='-', label='ALDH1A1')
plt.plot(df.index, df['ALDH1A2'], marker='s', linestyle='-', label='ALDH1A2')
plt.plot(df.index, df['ALDH1A3'], marker='^', linestyle='-', label='ALDH1A3')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Macula.tiff",
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
fig.set_size_inches(8,4)
plt.plot(df.index, df['ALDH1A1'], marker='o', linestyle='-', label='ALDH1A1')
plt.plot(df.index, df['ALDH1A2'], marker='s', linestyle='-', label='ALDH1A2')
plt.plot(df.index, df['ALDH1A3'], marker='^', linestyle='-', label='ALDH1A3')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Peripheral.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)


#########
markers = ['CYP26A1', 'CYP26B1','CYP26C1','TBX5','VAX2']
df = sc.get.obs_df(adata, markers + ["Region","Days"])
grouped_df = df.groupby(["Region","Days"])
result = grouped_df.mean().reset_index()
result.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege.csv")
df = result.loc[result.Region=="Macula",:]
# Convert the 'Date' column to datetime format

# Set the 'Date' column as the index
df.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8,4)
plt.plot(df.index, df['CYP26A1'], marker='o', linestyle='-', label='CYP26A1')
plt.plot(df.index, df['CYP26B1'], marker='s', linestyle='-', label='CYP26B1')
plt.plot(df.index, df['CYP26C1'], marker='^', linestyle='-', label='CYP26C1')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Macula2.tiff",
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
fig.set_size_inches(8,4)
plt.plot(df.index, df['CYP26A1'], marker='o', linestyle='-', label='CYP26A1')
plt.plot(df.index, df['CYP26B1'], marker='s', linestyle='-', label='CYP26B1')
plt.plot(df.index, df['CYP26C1'], marker='^', linestyle='-', label='CYP26C1')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Peripheral2.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

#########
markers = ['CYP26A1', 'CYP26B1','CYP26C1','TBX5','VAX2']
df = sc.get.obs_df(adata, markers + ["Region","Days"])
grouped_df = df.groupby(["Region","Days"])
result = grouped_df.mean().reset_index()
result.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege.csv")
df = result.loc[result.Region=="Macula",:]
# Convert the 'Date' column to datetime format

# Set the 'Date' column as the index
df.set_index('Days', inplace=True)

# Plot multiple lines using plt
plt.clf()
fig = plt.gcf()
fig.set_size_inches(8,4)
plt.plot(df.index, df['TBX5'], marker='s', linestyle='-', label='TBX5')
plt.plot(df.index, df['VAX2'], marker='^', linestyle='-', label='VAX2')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Macula3.tiff",
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
fig.set_size_inches(8,4)
plt.plot(df.index, df['TBX5'], marker='s', linestyle='-', label='TBX5')
plt.plot(df.index, df['VAX2'], marker='^', linestyle='-', label='VAX2')

# Customize the plot (optional)
plt.title('')
plt.xlabel('')
plt.ylabel('log CPM')
plt.grid(True)
plt.legend()  # Add legend to distinguish lines
plt.ylim(0, 2.5)
plt.tight_layout()
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/compare_pre_knowlege_Peripheral3.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)