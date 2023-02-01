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
df["scpred_prediction"] = df.majorclass.replace(
    {
        "AC Precursor": "AC",
        "BC Precursor": "BC",
        "Cone Precursor": "Cone",
        "GABAergic": "AC",
        "Glycinergic": "AC",
        "HC0": "HC",
        "HC1": "HC",
        "MG": "MG",
        "ML_Cone": "Cone",
        "NRPC": "RPC",
        "OFF-BC": "BC",
        "OFF_MGC": "RGC",
        "ON-BC": "BC",
        "ON_MGC": "RGC",
        "PRPC": "RPC",
        "RBC": "BC",
        "RGC Precursor": "RGC",
        "Rod": "Rod",
        "Rod Precursor": "Rod",
        "SACs": "AC",
        "S_Cone": "Cone",
        "dual ACs": "AC",
    }
)

df["Weeks"] = df.Days.map(
    {
        70: "FW10",
        79: "FW10",
        87: "FW13",
        91: "FW13",
        100: "FW16",
        103: "FW16",
        116: "FW16",
        120: "FW16",
        136: "FW20",
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
ax = counts.plot(x="Weeks", kind="barh", stacked=True, title=region)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.xlabel("")
plt.savefig(output_file_path + region + ".svg", bbox_inches="tight")