import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc

NRPC = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC.obs["subclass"] = np.nan

clusters = ["2", "13", "12"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["0", "1", "18"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["9", "10"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["4"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["17"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["7", "14"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "RGC"

df = NRPC.obs
df["majorclass"] = df["subclass"]
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
grouped = sample_df_Macula[["majorclass", "Days"]].groupby(["majorclass"])
mean = grouped.mean()
std = grouped.std()


def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Macula")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Macula.svg",transparent=True
    )


size = sample_df_Macula.groupby("majorclass").size() / np.sum(
    sample_df_Macula.groupby("majorclass").size()
)
plot_normal_distributions(mean.Days, std.Days, mean.index, size)

################################################################################################
sample_df_Peripheral = sample_df[sample_df.Region == "Peripheral"]
grouped = sample_df_Peripheral[["majorclass", "Days"]].groupby(["majorclass"])
mean = grouped.mean()
std = grouped.std()

size = sample_df_Peripheral.groupby("majorclass").size() / np.sum(
    sample_df_Peripheral.groupby("majorclass").size()
)


def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Peripheral")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Peripheral.svg",transparent=True
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)
