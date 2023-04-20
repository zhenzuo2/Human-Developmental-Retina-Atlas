import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

pd.options.display.max_columns = None
scv.set_figure_params(dpi=600, dpi_save=600)

adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)
adata_result

Time = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
)
adata_result = adata_result[Time["Unnamed: 0"]]

adata_result.obs["latent_time"] = Time["latent_time"].values
modules = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/hm.get_yticklabels.csv"
)
sc.pp.normalize_total(adata_result, target_sum=1e4, exclude_highly_expressed=True)
sc.pp.log1p(adata_result)

res = []
for x in [1, 2, 3, 4]:
    temp = []
    for gene in modules[modules.Module == x]["Unnamed: 0"][:100]:
        gs = [x[0] for x in adata_result[:, gene].X.toarray()]
        gs = gs / (np.sum(gs))
        time = adata_result.obs["latent_time"]
        temp = temp + [np.sum(gs * time)]
    res = res + [temp]

# create a box plot with different colors for each list
ax = sns.boxplot(
    data=[res[i] for i in [2, 1, 3, 0]],
    palette=[
        "#1F77B4",
        "#FF7F0E",
        "#2CA02C",
        "#D62728",
        "#9467BD",
        "#8C564B",
        "#E377C2",
        "#7F7F7F",
        "#BCBD22",
        "#17BECF",
    ],
)

# set the title and axis labels
# set the title and y-axis label
ax.set_title("Latent time of Gene Expression Groupped by Module")
ax.set_ylabel("Gene Expression Weighted Latent Time")

ax.set_xticklabels(
    [["Module 1", "Module 2", "Module 3", "Module 4"][i] for i in [2, 1, 3, 0]]
)

# show the plot
plt.savefig(
    "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/Boxplot_module.svg",
    bbox_inches="tight",
)
plt.clf()

n = 20
adata_result.obs["latent_time_grouped"] = pd.qcut(adata_result.obs.latent_time, q=n)

for module_index in [1, 2, 3, 4]:
    matrix = adata_result[
        :, modules[modules.Module == module_index]["Unnamed: 0"][:50]
    ].X.toarray()

    # Create a list of row group labels for the matrix
    labels = adata_result.obs["latent_time_grouped"].values

    # Convert the matrix to a pandas dataframe and add the labels as a new column
    df = pd.DataFrame(matrix)
    df["Group"] = labels

    # Group the dataframe by the 'Group' column and calculate the mean for each group
    grouped = df.groupby("Group").mean()

    # Remove the 'Group' column from the grouped dataframe
    # grouped = grouped.drop('Group', axis=1)

    # Convert the grouped dataframe back to a numpy array
    grouped_matrix = grouped.to_numpy()

    # Create a 5x10 matrix with random values between 0 and 1
    matrix = grouped_matrix

    # Create an array of x values from 1 to 10
    x = np.arange(1, n + 1)

    # Loop through each column in the matrix and plot a line for it
    for i in range(matrix.shape[1]):
        y = matrix[:, i]
        y = stats.zscore(y)
        plt.plot(x, y)

    # Add labels and a title to the plot
    plt.xlabel("Latent Time (Bins)")
    plt.ylabel("Gene Expression Z score")
    plt.title("Module " + str(module_index))

    # Show the plot
    plt.savefig(
        "/storage/singlecell/zz4/fetal_bash/results/PRPC_hotspot/Module_"
        + str(module_index)
        + "_gene_expression_curve.svg"
    )
    plt.clf()

adata_result_Macula = adata_result[adata_result.obs.Region == "Macula"]
adata_result_Peripheral = adata_result[adata_result.obs.Region == "Peripheral"]
res_Macula = []
for x in [1, 2, 3, 4]:
    temp = []
    for gene in modules[modules.Module == x]["Unnamed: 0"]:
        gs = [x[0] for x in adata_result_Macula[:, gene].X.toarray()]
        gs = gs / (np.sum(gs))
        time = adata_result_Macula.obs["Days"]
        temp = temp + [np.sum(gs * time)]
    res_Macula = res_Macula + [temp]
res_Peripheral = []
for x in [1, 2, 3, 4]:
    temp = []
    for gene in modules[modules.Module == x]["Unnamed: 0"][:50]:
        gs = [x[0] for x in adata_result_Peripheral[:, gene].X.toarray()]
        gs = gs / (np.sum(gs))
        time = adata_result_Peripheral.obs["Days"]
        temp = temp + [np.sum(gs * time)]
    res_Peripheral = res_Peripheral + [temp]
data = pd.DataFrame(
    {
        "Days": res_Macula[0]
        + res_Macula[1]
        + res_Macula[2]
        + res_Macula[3]
        + res_Peripheral[0]
        + res_Peripheral[1]
        + res_Peripheral[2]
        + res_Peripheral[3],
        "Module": ["Module1"] * len(res_Macula[0])
        + ["Module2"] * len(res_Macula[1])
        + ["Module3"] * len(res_Macula[2])
        + ["Module4"] * len(res_Macula[3])
        + ["Module1"] * len(res_Peripheral[0])
        + ["Module2"] * len(res_Peripheral[1])
        + ["Module3"] * len(res_Peripheral[2])
        + ["Module4"] * len(res_Peripheral[3]),
        "Region": ["Macula"] * len(res_Macula[0])
        + ["Macula"] * len(res_Macula[1])
        + ["Macula"] * len(res_Macula[2])
        + ["Macula"] * len(res_Macula[3])
        + ["Peripheral"] * len(res_Peripheral[0])
        + ["Peripheral"] * len(res_Peripheral[1])
        + ["Peripheral"] * len(res_Peripheral[2])
        + ["Peripheral"] * len(res_Peripheral[3]),
    }
)
sns.boxplot(
    x="Module",
    y="Days",
    hue="Region",
    data=data,
    palette=["#EF553B", "#636EFA"],
    order=["Module3", "Module2", "Module4", "Module1"],
    showfliers=False,
)
sns.despine(offset=10, trim=True)
plt.savefig("sns.boxplot.svg")