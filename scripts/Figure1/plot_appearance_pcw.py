import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import scanpy as sc
import numpy as np

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata = adata[adata.obs.Region.isin(["Macula", "Peripheral"])]
adata.obs["Region"] = adata.obs["Region"].replace({"Peripheral": "Periphery"})
df = adata.obs
sample_size = 1000

# group the dataframe by the 'Group' column
groups = df.groupby("Days")

# create an empty dataframe to store the sample data
sample_df = pd.DataFrame(columns=df.columns)

# sample equal number of rows by group with replicates
for group, data in groups:
    sample_data = data.sample(n=sample_size, replace=True)
    sample_df = pd.concat([sample_df, sample_data])

sample_df.majorclass = sample_df.majorclass.astype(str)

sample_df


def confidence_interval_length(std, n):
    # Degrees of freedom
    df = n - 1
    # t-value for a 95% confidence interval (two-tailed)
    t_value = 1.96  # approximate value for large sample sizes
    # Calculate standard error
    stderr = std / np.sqrt(n)
    # Calculate length of confidence interval (one side)
    interval_length = t_value * stderr
    return interval_length


# Group by 'Region' and 'majorclass', calculating mean and std
group_stats = sample_df.groupby(["Region", "majorclass"]).agg(
    {"Days": ["mean", "std", "count"]}
)
group_stats["Days", "CI_length"] = confidence_interval_length(
    group_stats["Days", "std"], group_stats["Days", "count"]
)
group_stats = group_stats.reset_index()
group_stats.columns = ["Region", "majorclass", "mean", "std", "count", "CI_length"]

group_stats.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/plot_appearance_pcw.csv"
)
colors = {"Macula": "#F8766D", "Periphery": "#00BFC4"}
# Iterate over unique regions
for region in group_stats["Region"].unique():
    region_data = group_stats[group_stats["Region"] == region]
    region_data.index = region_data.majorclass.values
    region_data = region_data.loc[["RGC", "HC", "Cone", "AC", "Rod", "BC", "MG"],]
    plt.plot(
        region_data["majorclass"],
        region_data["mean"],
        marker="o",
        label=region,
        color=colors[region],
    )

    plt.errorbar(
        region_data["majorclass"],
        region_data["mean"],
        yerr=region_data["CI_length"],  # Use the std column for error bars
        # marker="o",
        fmt="o",
        # label=region,
        ecolor=colors[region],
        color=colors[region],
        capsize=10,
    )

#plt.title("Mean Days after Conception by Major Class")
plt.title("")  
plt.xlabel("Major Class")
plt.ylabel("Mean Days after Conception")
plt.xticks(rotation=45)
plt.legend()
plt.grid(True)
plt.tight_layout()
fig = plt.gcf()
fig.set_size_inches(3, 3)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/plot_appearance_pcw.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi=600,
)
