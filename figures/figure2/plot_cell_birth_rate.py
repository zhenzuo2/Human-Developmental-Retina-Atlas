import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
NRPC = sc.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
NRPC = NRPC[NRPC.obs.leiden!="12"]
NRPC.obs["subclass"] = np.nan

sc.pp.neighbors(NRPC, use_rep="X_scVI")
sc.tl.leiden(NRPC, resolution=3)

clusters = ["28", "3", "13", "15"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "AC"
clusters = ["6", "30","35","25"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "HC"
clusters = ["21","14","27"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "BC"
clusters = ["33","0"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Rod"
clusters = ["12","22"]
NRPC.obs.loc[NRPC.obs.leiden.isin(clusters), "subclass"] = "Cone"
clusters = ["4", "19","23"]
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
    plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Macula.svg",transparent=True,bbox_inches='tight'
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
    plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.savefig(
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Peripheral.svg",transparent=True,bbox_inches='tight'
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)
