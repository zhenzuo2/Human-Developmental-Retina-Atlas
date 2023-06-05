import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/storage/singlecell/zz4/fetal_snakemake/results/pseudotime/PRPC_MG_latent_time.csv"
)
output_file_path = "/storage/singlecell/zz4/fetal_snakemake/figures/figure6/"
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

sns.boxplot(
    x="Weeks",
    y="latent_time",
    hue="Region",
    data=df,
    order=["PCW10", "PCW13", "PCW16", "PCW20", "PCW23"],
    showfliers=False,
)
plt.ylabel("Gene-shared latent time")
plt.xlabel("Post Conception Week")
plt.legend(
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.savefig(
    output_file_path + "latent_time_WEEK.svg",
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

sns.boxplot(x="Days", y="latent_time", hue="Region", data=df, showfliers=False)
plt.ylabel("Gene-shared latent time")
plt.xlabel("Post Conception Days")
plt.legend(
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.savefig(
    output_file_path + "latent_time_Days.svg",
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)

sns.boxplot(x="majorclass", y="latent_time", data=df, showfliers=False)
plt.ylabel("Differentiation Potency")
plt.xlabel("Cell Type")
plt.legend(
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    borderaxespad=0,
)
plt.savefig(
    output_file_path + "latent_time_Region.svg",
    bbox_inches="tight",
    transparent=True,
    backend="cairo",
)
