import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({"font.size": 18})
NRPC = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.h5ad"
)
NRPC = NRPC[NRPC.obs.subclass != "NRPC"]
df = NRPC.obs
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
cols = {
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "AC": "#8c564b",
}

def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i], color=cols[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Macula")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    ax.set_ylim([0, 0.006])
    ax.set_xlim([70, 160])
    for line in leg.get_lines():
        line.set_linewidth(10.0)
    plt.legend('',frameon=False)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate in the Macula.tiff",
        transparent=True,
        bbox_inches="tight",
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

mean.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate mean.csv")
std.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate std.csv")


size = sample_df_Peripheral.groupby("subclass").size() / np.sum(
    sample_df_Peripheral.groupby("subclass").size()
)
size.to_csv("/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate size.csv")

def plot_normal_distributions(means, stds, labels, size):
    x = np.linspace(50, 180, 1000)
    i = 0
    fig, ax = plt.subplots()
    for mean, std in zip(means, stds):
        y = (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        y = y * size[i]
        ax.plot(x, y, label=labels[i], color=cols[labels[i]])
        i = i + 1
    # get the legend object
    leg = ax.legend()
    ax.set_title("Cell Birth Rate in the Peripheral")
    ax.set_xlabel("Post Conception Days")
    ax.set_ylabel("Normalized Birth Rate")
    ax.set_ylim([0, 0.006])
    ax.set_xlim([70, 160])
    for line in leg.get_lines():
        line.set_linewidth(10.0)
    plt.legend('',frameon=False)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate in the Peripheral.tiff",
        transparent=True,
        bbox_inches="tight",
    )


plot_normal_distributions(mean.Days, std.Days, mean.index, size)

# Define your order and colors
order = ['RGC', 'Cone', 'HC', 'AC', 'BC', 'Rod']
cols = {
    "Rod": "#17becf",
    "BC": "#bcbd22",
    "RGC": "#d62728",
    "Cone": "#e377c2",
    "HC": "#2ca02c",
    "AC": "#8c564b",
}

# Create custom legend handles in the specified order
legend_handles = [Line2D([0], [0], color=cols[label], lw=8, label=label) for label in order]

# Plot the legend with a single column (horizontal)
fig, ax = plt.subplots()
ax.legend(handles=legend_handles, loc='center', bbox_to_anchor=(0.5, -0.2), ncol=1)

# Hide the axes
ax.set_axis_off()

# Show the plot
plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/NRPC Birth Rate legend.tiff",
        transparent=True,
        bbox_inches="tight",
    )
