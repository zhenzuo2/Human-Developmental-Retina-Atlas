import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import matplotlib.ticker as mtick

input_file_path = "/storage/singlecell/zz4/fetal_snakemake/results/cell_annotation_results/filtered_major_class_MG.csv"
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/"
region = "Macula"

try:
    os.makedirs(output_file_path)
except FileExistsError:
    pass

df = pd.read_csv(input_file_path)
df = df.loc[df.Region == region,]
df = df.reset_index()

df["Weeks"] = df.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)

df["majorclass"] = df.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)

counts = pd.DataFrame(np.zeros([5, 8]))
counts.columns = [
    "RPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]
counts.index = ["PCW10", "PCW13", "PCW16", "PCW20", "PCW23"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.majorclass[i]] = (
        counts.at[df.Weeks[i], df.majorclass[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
# ax.spines['bottom'].set_visible(False)
ax.spines["left"].set_visible(False)

# ax.get_xaxis().set_ticks([])
# ax.get_yaxis().set_ticks([])

plt.legend(
    bbox_to_anchor=(1, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.xlabel("")
plt.ylabel("Post Conception Week")
plt.savefig(output_file_path + region + ".svg", bbox_inches="tight")

############################################################################################################################################
region = "Peripheral"

try:
    os.makedirs(output_file_path)
except FileExistsError:
    pass

df = pd.read_csv(input_file_path)
df = df.loc[df.Region == region,]
df = df.reset_index()

df["Weeks"] = df.Days.map(
    {
        70: "PCW10",
        79: "PCW10",
        87: "PCW13",
        91: "PCW13",
        100: "PCW16",
        103: "PCW16",
        116: "PCW16",
        120: "PCW16",
        136: "PCW20",
        137: "PCW20",
        141: "PCW20",
        142: "PCW20",
        162: "PCW23",
        165: "PCW23",
    }
)

df["majorclass"] = df.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)

counts = pd.DataFrame(np.zeros([5, 8]))
counts.columns = [
    "RPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]
counts.index = ["PCW10", "PCW13", "PCW16", "PCW20", "PCW23"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.majorclass[i]] = (
        counts.at[df.Weeks[i], df.majorclass[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
# ax.spines['bottom'].set_visible(False)
ax.spines["left"].set_visible(False)

# ax.get_xaxis().set_ticks([])
# ax.get_yaxis().set_ticks([])

plt.legend(
    bbox_to_anchor=(1, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.xlabel("")
plt.ylabel("Post Conception Week")
plt.savefig(output_file_path + region + ".svg", bbox_inches="tight")
