import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scvelo as scv

palette={
    "NRPC": "#9467bd",
    "AC Precursor": "#17becf",
    "GABAergic": "#bcbd22",
    "Glycinergic": "#d62728",
    "SACs": "#ff7f0e",
    "dual ACs": "#e377c2",
}
adata = sc.read("/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad")

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

df = adata.obs
df = df[df.Days > 0]

# define the sample size
sample_size = 20000

# group the dataframe by the 'Group' column
groups = df.groupby("Days")

# create an empty dataframe to store the sample data
sample_df = pd.DataFrame(columns=df.columns)

# sample equal number of rows by group with replicates
for group, data in groups:
    sample_data = data.sample(n=sample_size, replace=True)
    sample_df = pd.concat([sample_df, sample_data])

sample_df_Macula = sample_df[sample_df.Region == "Macula"]
grouped = sample_df_Macula[["subclass", "Days"]].groupby(["subclass"])
mean = grouped.mean()
std = grouped.std()


def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i],color = palette[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Macula")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/Cell Birth Rate in the Macula_AC.svg",transparent=True
    )


size = sample_df_Macula.groupby("subclass").size() / np.sum(
    sample_df_Macula.groupby("subclass").size()
)
plot_normal_distributions(mean.Days, std.Days, mean.index, size)

################################################################################################
sample_df_Peripheral = sample_df[sample_df.Region == "Peripheral"]
grouped = sample_df_Peripheral[["subclass", "Days"]].groupby(["subclass"])
mean = grouped.mean()
std = grouped.std()

size = sample_df_Peripheral.groupby("subclass").size() / np.sum(
    sample_df_Peripheral.groupby("subclass").size()
)


def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i],color = palette[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Peripheral")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/Cell Birth Rate in the Peripheral_AC.svg",transparent=True
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)