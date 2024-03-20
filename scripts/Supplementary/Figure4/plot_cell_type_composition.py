import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import matplotlib.ticker as mtick

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

input_file_path = "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv"
output_file_path = "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure4/"
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

counts = pd.DataFrame(np.zeros([5, 9]))
counts.columns = [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]
counts.index = ["PCW10", "PCW13", "PCW15", "PCW19", "PCW23"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.majorclass[i]] = (
        counts.at[df.Weeks[i], df.majorclass[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
counts.to_csv(output_file_path + region + "_cell_type_composition.csv")
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region,color = [pairwise_dict[x] for x in [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]])
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
plt.savefig(output_file_path + region + "_cell_type_composition.svg", bbox_inches="tight")

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

counts = pd.DataFrame(np.zeros([5, 9]))
counts.columns = [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]
counts.index = ["PCW10", "PCW13", "PCW15", "PCW19", "PCW23"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.majorclass[i]] = (
        counts.at[df.Weeks[i], df.majorclass[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
counts.to_csv(output_file_path + region + "_cell_type_composition.csv")
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region,color = [pairwise_dict[x] for x in [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]])
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
plt.savefig(output_file_path + region + "_cell_type_composition.svg", bbox_inches="tight")

############################################################################################################################################
region = "Whole Eye"

try:
    os.makedirs(output_file_path)
except FileExistsError:
    pass

df = pd.read_csv(input_file_path)
df = df.loc[df.Region == region,]
df = df.reset_index()

df["Weeks"] = df.Days.map(
    {
        59: "PCW8",
    }
)

counts = pd.DataFrame(np.zeros([1, 9]))
counts.columns = [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]
counts.index = ["PCW8"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.majorclass[i]] = (
        counts.at[df.Weeks[i], df.majorclass[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
counts.to_csv(output_file_path + region + "_cell_type_composition.csv")
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region,color = [pairwise_dict[x] for x in [
    "PRPC",
    "NRPC",
    "RGC",
    "Cone",
    "HC",
    "AC",
    "Rod",
    "BC",
    "MG",
]])
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
plt.savefig(output_file_path + region + "_cell_type_composition.svg", bbox_inches="tight")

