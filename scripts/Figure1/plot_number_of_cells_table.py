import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import table
from PIL import Image
import matplotlib
matplotlib.font_manager._load_fontmanager(try_read_cache=False)

plt.rcParams["font.family"] = "Arial"
plt.rcParams['font.size'] = 20

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
adata.obs["Weeks"] = adata.obs.Days.map(
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

df = pd.DataFrame(adata.obs["Weeks"].value_counts()).loc[
    ["PCW8", "PCW10", "PCW13", "PCW15", "PCW19", "PCW23"], :
]
df["Location"] = [
    "Whole Eye",
    "Macula and Periphery",
    "Macula and Periphery",
    "Macula and Periphery",
    "Macula and Periphery",
    "Macula and Periphery",
]
df["Method"] = "Sn-Multiomics"
df.columns = ["Number of Nuclei After QC", "Location", "Method"]
df["Number of Nuclei After QC"] = df.apply(lambda x: "{:,}".format(x["Number of Nuclei After QC"]), axis=1)
df = df.loc[:,["Location", "Method","Number of Nuclei After QC"]]

def draw_table_with_colors(df, output_path):
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(5, 4))  # Adjust the figure size as needed
    # Plot an empty table
    ax.axis("off")
    tab = table(
        ax,
        df,
        loc="center",
        cellLoc="center",
        cellColours=[["#ffffff"] * len(df.columns)] * len(df.index),
        colWidths=list([0.5] * len(df.index)),
    )
    # Set header color
    for key, cell in tab.get_celld().items():
        if key[0] == 0:
            cell.set_text_props(weight="bold", color="white",size = 30)
            cell.set_facecolor("#3498db")  # Blue color for header
    #plt.tight_layout()
    # Save the figure as a TIFF image
    plt.savefig(
        output_path,
        format="tiff",
        bbox_inches="tight",
        pad_inches=0.5,
        dpi=300,
        transparent=True,
    )
    # Close the figure to free up resources
    plt.close(fig)


draw_table_with_colors(
    df,
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/Numer of cells each PCW.tiff",
)
