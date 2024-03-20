# Import packages
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib
from PIL import Image

matplotlib.font_manager._load_fontmanager(try_read_cache=False)
plt.rcParams["font.family"] = "Arial"

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/AC_w_NRPC.h5ad"
)

adata.obs["subclass"] = pd.Categorical(
    list(adata.obs["subclass"]),
    categories=["NRPC", "AC Precursor", "GABAergic", "Glycinergic", "SACs", "dual ACs"],
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

cells = []
for sp in set(adata.obs.subclass):
    cells = cells + [adata[adata.obs.subclass == sp].obs.index[0]]

for region in set(adata.obs.Region):
    for weeks in ['PCW10', 'PCW13', 'PCW15', 'PCW19', 'PCW23',]:
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "subclass"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "subclass"]

        sc.pl.umap(
            adata,
            color="temp",
            size=40,
            title="",
            frameon=False,
            legend_loc=None,
            palette={
                "NRPC": "#9467bd",
                "AC Precursor": "#17becf",
                "GABAergic": "#bcbd22",
                "Glycinergic": "#d62728",
                "SACs": "#ff7f0e",
                "dual ACs": "#e377c2",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.savefig(
            "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_"
            + region
            + "_"
            + weeks
            + ".png",
            dpi=300,
            bbox_inches="tight",
            transparent=True,
        )

region = "Whole Eye"
weeks = "PCW8"
adata.obs["temp"] = np.nan
subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "subclass"]
adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "subclass"]
sc.pl.umap(
    adata,
    color="temp",
    size=40,
    title="",
    frameon=False,
    legend_loc=None,
    palette={
        "NRPC": "#9467bd",
        "AC Precursor": "#17becf",
        "GABAergic": "#bcbd22",
        "Glycinergic": "#d62728",
        "SACs": "#ff7f0e",
        "dual ACs": "#e377c2",
    },
)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_"
    + 'Whole Eye'
    + "_"
    + "PCW8"
    + ".png",
    dpi=300,
    bbox_inches="tight",
    transparent=True,
)

def merge_images(file_paths):
    # Initialize an empty list to store image objects
    images = []

    # Loop through each sublist in the double list
    for sublist in file_paths:
        # Initialize an empty list to store image arrays
        image_arrays = []

        # Loop through each file path in the sublist
        for file_path in sublist:
            # Open the image using PIL
            img = Image.open(file_path)

            # Convert the image to an array
            img_array = np.array(img)

            # Append the image array to the list
            image_arrays.append(img_array)

        # Stack the image arrays horizontally to create a row of images
        row_images = np.hstack(image_arrays)

        # Append the row of images to the list
        images.append(row_images)

    # Stack the rows of images vertically to create the final merged image
    merged_image = np.vstack(images)

    # Display the merged image using Matplotlib
    plt.imshow(merged_image)
    plt.axis("off")
    plt.show()

    return merged_image


# Example usage
file_paths = [
    [
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Peripheral_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Peripheral_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Peripheral_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Peripheral_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Peripheral_PCW23.png",
    ],
    [
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Macula_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Macula_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Macula_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Macula_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_Macula_PCW23.png",
    ],
]

merged_image = merge_images(file_paths)
im = Image.fromarray(merged_image)
im.save(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/overall_umap_by_subclass_PCW.png",
    dpi=(300, 300),
)
#######################################################################################################################

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Define the color dictionary
color_dict = {
    "NRPC": "#9467bd",
    "AC Precursor": "#17becf",
    "GABAergic": "#bcbd22",
    "Glycinergic": "#d62728",
    "SACs": "#ff7f0e",
    "dual ACs": "#e377c2",
}
# Set the desired marker size
markersize = 15

# Create custom legend handles with solid dots of specified size
legend_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=markersize) for color in color_dict.values()]
legend_labels = list(color_dict.keys())

# Create the legend
plt.legend(legend_handles, legend_labels, loc='center', bbox_to_anchor=(0.5, -0.2), ncol=len(color_dict), frameon=False)

# Remove axis ticks and labels
plt.xticks([])
plt.yticks([])

# Save the legend as a TIFF file
plt.savefig('/storage/chentemp/zz4/adult_dev_compare/figures/Figure4/legend_bigger_dots.tiff', bbox_inches='tight', pad_inches=0.1, format='tiff',dpi = 600)

# Show the plot
plt.show()