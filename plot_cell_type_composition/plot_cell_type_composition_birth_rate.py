import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import matplotlib.ticker as mtick

input_file_path = "/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/merged_raw_filtered_umap_10000_major_sub_class.obs.csv"
region = "Macula"

df = pd.read_csv(input_file_path)
df = df.loc[df.Region == region,]
df = df.loc[
    df.Days.isin(
        [
            70.0,
            79.0,
            87.0,
            91.0,
            100.0,
            103.0,
            # 116.0,
            120.0,
            136.0,
            # 137.0,
            141.0,
            142.0,
            162.0,
            165.0,
        ]
    )
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
counts = pd.DataFrame(np.zeros([12, 8]))
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
counts.index = [
    70.0,
    79.0,
    87.0,
    91.0,
    100.0,
    103.0,
    # 116.0,
    120.0,
    136.0,
    # 137.0,
    141.0,
    142.0,
    162.0,
    165.0,
]
for i in range(df.shape[0]):
    counts.at[df.Days[i], df.scpred_prediction[i]] = (
        counts.at[df.Days[i], df.scpred_prediction[i]] + 1
    )
counts = counts.div(counts.sum(axis=1), axis=0) * 100
counts["Days"] = counts.index.values
ax = counts.plot(x="Days", kind="barh", stacked=True, title=region)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.xlabel("")
plt.savefig(region + ".svg", bbox_inches="tight")

import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import statsmodels.api as sm

# create sample dataframe
df = counts.iloc[:, :8].diff().iloc[1:, :]

# plot a line chart with regression curves
fig, ax = plt.subplots()
for col in df.columns:
    x_data = df.index
    y_data = df[col]
    yhat = sm.nonparametric.lowess(y_data, x_data, return_sorted=False)
    ax.plot(x_data, yhat, label=col)
ax.legend()

# display the plot
plt.legend(loc=(1.04, 0))
plt.title('Cell Type Proportions Changes From Last Time Point (Macula)')
plt.xlabel('Days')
plt.ylabel('Cell Type Proportions Changes')
plt.savefig(region + "2.svg",bbox_inches='tight')