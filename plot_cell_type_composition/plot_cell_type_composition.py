import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import matplotlib.ticker as mtick

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
region = sys.argv[3]

try:
    os.makedirs(output_file_path)
except FileExistsError:
    pass

df = pd.read_csv(input_file_path)
df = df.loc[
    df.Region == region,
]
df = df.reset_index()
df["Weeks"] = df.Days.map(
    {
        70: "FW10",
        79: "FW10",
        87: "FW13",
        91: "FW13",
        100: "FW16",
        103: "FW16",
        116: "FW16",
        137: "FW20",
        141: "FW20",
        142: "FW20",
        162: "FW23",
        165: "FW23",
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
counts.index = ["FW10", "FW13", "FW16", "FW20", "FW23"]
for i in range(df.shape[0]):
    counts.at[df.Weeks[i], df.scpred_prediction[i]] = (
        counts.at[df.Weeks[i], df.scpred_prediction[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Weeks"] = counts.index.values
ax = counts.plot(x="Weeks", kind="bar", stacked=True, title=region)
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.xlabel("")
plt.savefig(output_file_path + region + ".svg", bbox_inches='tight')