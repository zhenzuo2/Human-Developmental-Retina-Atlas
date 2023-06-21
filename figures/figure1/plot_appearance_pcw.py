import scvelo as scv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pandas as pd
adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
df = adata.obs
sample_size = 20000

# group the dataframe by the 'Group' column
groups = df.groupby("Days")

# create an empty dataframe to store the sample data
sample_df = pd.DataFrame(columns=df.columns)

# sample equal number of rows by group with replicates
for group, data in groups:
    sample_data = data.sample(n=sample_size, replace=True)
    sample_df = pd.concat([sample_df, sample_data])
sample_df.majorclass = sample_df.majorclass.astype(str)

cell_types = ["RGC","HC","AC","Cone","Rod","BC","MG"]
for cell_type in cell_types:
    df = sample_df[sample_df.majorclass==cell_type].groupby(["Region","majorclass","Days"])["sampleid"].count()
    total = sample_df[sample_df.majorclass==cell_type].groupby(["majorclass"])["sampleid"].count().values[0]
    n=int(len(df)/2.0)
    def find_index_10_quantile(numbers,total_sum):
        #total_sum = sum(numbers)
        quantile_10 = total_sum * 0.05
        aggregated_sum = 0

        for i, num in enumerate(numbers):
            aggregated_sum += num
            if aggregated_sum >= quantile_10:
                return i

        return -1  # If the quantile is not reached, return -1
    # Example usage:
    numbers = df.values[:n]
    index_10_quantile = find_index_10_quantile(numbers,total)
    print(cell_type)
    print("Index at which the 10% quantile is reached:", [x[2] for x in df.index][index_10_quantile])
    # Example usage:
    numbers = df.values[n:]
    index_10_quantile = find_index_10_quantile(numbers,total)
    print("Index at which the 10% quantile is reached:", [x[2] for x in df.index][index_10_quantile])

RGC = [70,70]
HC = [70,79]
AC = [91,100]
Cone = [91,100]
Rod = [103,116]
BC = [103,142]
MG = [116,165]

matplotlib.rcParams.update({'font.size': 22})
# create data
x = [x/7 for x in [70,70,91,91,103,103,116]]
y = [x/7 for x in [70,79,100,100,116,142,165]]
  
# plot lines
plt.plot(["RGC","HC","Cone","AC","Rod","BC","MG"], x, label = "Macula", marker=".", markersize=20,color = "#F8766D")
plt.plot(["RGC","HC","Cone","AC","Rod","BC","MG"], y, label = "Peripheral", marker=".", markersize=20,color = "#00BFC4")
plt.xlabel("Cell Type")
plt.ylabel("Appearance PWC")
plt.plot(["RGC","RGC"], [x/7 for x in RGC], linestyle='dashed',color = "black")
plt.plot(["AC","AC"], [x/7 for x in AC], linestyle='dashed',color = "black")
plt.plot(["HC","HC"], [x/7 for x in HC], linestyle='dashed',color = "black")
plt.plot(["Cone","Cone"], [x/7 for x in Cone], linestyle='dashed',color = "black")
plt.plot(["Rod","Rod"], [x/7 for x in Rod], linestyle='dashed',color = "black")
plt.plot(["BC","BC"], [x/7 for x in BC], linestyle='dashed',color = "black")
plt.plot(["MG","MG"], [x/7 for x in MG], linestyle='dashed',color = "black")
plt.legend()
plt.xticks(rotation=90)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure1/plot_appearance_pcw.svg",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
    backend="cairo",
)