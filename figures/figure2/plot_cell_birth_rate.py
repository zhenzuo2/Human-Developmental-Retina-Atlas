import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
)
mapping = {"NRPC": "RPC", "PRPC": "RPC"}
df["majorclass"] = df["majorclass"].replace(mapping)

df["majorclass"] = pd.Categorical(
    list(df["majorclass"]),
    categories=[
        "RPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
        "MG",
    ],
)
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
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Macula.svg"
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
        "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/Cell Birth Rate in the Peripheral.svg"
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)
