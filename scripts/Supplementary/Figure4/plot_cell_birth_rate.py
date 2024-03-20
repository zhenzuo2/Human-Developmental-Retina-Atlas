import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
pairwise_dict = {
    "MG": "#9467bd",
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "NRPC": "#ff7f0e",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "PRPC": "#1f77b4",
    "AC": "#8c564b",
}

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.majorclass.isin(['AC', 'BC', 'Cone', 'HC', 'MG', 'RGC', 'Rod'])]
df = adata.obs

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
        ax.plot(x, y, label=labels[i],color = pairwise_dict[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Macula")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    ax.set_ylim([0, 0.012])
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Macula.svg",transparent=True,bbox_inches='tight'
    )


size = sample_df_Macula.groupby("majorclass").size() / np.sum(
    sample_df_Macula.groupby("majorclass").size()
)

mean.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Macula mean.csv")
std.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Macula std.csv")
size.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Macula size.csv")
plot_normal_distributions(mean.Days, std.Days, mean.index, size)

################################################################################################
sample_df_Peripheral = sample_df[sample_df.Region == "Peripheral"]
grouped = sample_df_Peripheral[["majorclass", "Days"]].groupby(["majorclass"])
mean = grouped.mean()
std = grouped.std()

size = sample_df_Peripheral.groupby("majorclass").size() / np.sum(
    sample_df_Peripheral.groupby("majorclass").size()
)

mean.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Peripheral mean.csv")
std.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Peripheral std.csv")
size.to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Peripheral size.csv")

def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i],color = pairwise_dict[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Peripheral")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    ax.set_ylim([0, 0.012])
    for line in leg.get_lines():
        line.set_linewidth(5.0)
    plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/Cell Birth Rate in the Peripheral.svg",transparent=True,bbox_inches='tight'
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)