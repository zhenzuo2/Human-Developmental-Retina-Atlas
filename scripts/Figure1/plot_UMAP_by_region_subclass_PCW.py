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

color_panels = [
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
    "#b15928",
    "#a6cee3",
    "#b2df8a",
    "#fb9a99",
    "#fdbf6f",
    "#cab2d6",
    "#ffff99",
    "#b4c7e7",
    "#c7e9c0",
    "#fdc086",
    "#d6604d",
    "#9e9ac8",
    "#fdae61",
    "#5e3c99",
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#666666",
    "#1a1a1a",
    "#636363",
    "#bdbdbd",
]


adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

adata.obs["majorclass"] = pd.Categorical(
    list(adata.obs["majorclass"]),
    categories=[
        "PRPC",
        "NRPC",
        "RGC",
        "Cone",
        "HC",
        "AC",
        "Rod",
        "BC",
        "MG",
    ],
)
#############################################################################################################################
sc.pl.umap(
    adata,
    color="Region",
    size=5,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_region.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

adata.obs["Days"] = adata.obs["Days"].astype(float)
sc.pl.umap(
    adata,
    color="Days",
    size=5,
    legend_loc="right margin",
    title="",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_days.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.rcParams.update({"font.size": 12})
sc.pl.umap(
    adata,
    color="subclass",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
    palette=color_panels,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_subclass_legend.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

sc.pl.umap(
    adata,
    color="subclass",
    size=5,
    legend_loc=None,
    title="",
    return_fig=True,
    frameon=False,
    legend_fontoutline=True,
    palette=color_panels,
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_subclass.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
)

plt.rcParams.update({"font.size": 20})
sc.pl.umap(
    adata,
    color="majorclass",
    size=5,
    legend_loc="on data",
    title="",
    return_fig=True,
    frameon=False,
    legend_fontweight="bold",
    legend_fontoutline=True,
    palette="tab20",
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass.tiff",
    dpi=600,
    bbox_inches="tight",
    transparent=True,
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
for sp in set(adata.obs.majorclass):
    cells = cells + [adata[adata.obs.majorclass == sp].obs.index[0]]

for region in set(adata.obs.Region):
    for weeks in set(adata.obs.Weeks):
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "majorclass"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "majorclass"]
        adata.obs["temp"] = pd.Categorical(
            list(adata.obs["temp"]),
            categories=[
                "PRPC",
                "NRPC",
                "RGC",
                "Cone",
                "HC",
                "AC",
                "Rod",
                "BC",
                "MG",
            ],
        )
        sc.pl.umap(
            adata,
            color="temp",
            size=2,
            title="",
            frameon=False,
            legend_loc="None",
            palette={
                "MG": "#9467bd",
                "Rod": "#17becf",
                "BC": "#bcbd22",
                "RGC": "#d62728",
                "NRPC": "#ff7f0e",
                "Cone": "#e377c2",
                "HC": "#2ca02c",
                "PRPC": "#1f77b4",
                "AC": "#8c564b",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(5, 5)
        plt.savefig(
            "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_"
            + region
            + "_"
            + weeks
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
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Peripheral_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Peripheral_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Peripheral_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Peripheral_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Peripheral_PCW23.png",
    ],
    [
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Macula_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Macula_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Macula_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Macula_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_Macula_PCW23.png",
    ],
]

merged_image = merge_images(file_paths)
im = Image.fromarray(merged_image)
im.save(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure1/overall_umap_by_majorclass_PCW.png",
    dpi=(300, 300),
)
