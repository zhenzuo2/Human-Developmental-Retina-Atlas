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

adata = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/NRPC_annotated_terminal_states.h5ad"
)
adata_result = sc.read_h5ad(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/NRPC_res2.h5ad"
)

adata.obsm["X_umap"] = adata_result[adata.obs.index].obsm["X_umap"]

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
for sp in ['AC', 'HC', 'Cone', 'Rod', 'RGC', 'BC']:
    cells = cells + [adata[adata.obs.subclass == sp].obs.index[0]]

for region in set(adata.obs.Region):
    for weeks in set(adata.obs.Weeks):
        adata.obs["temp"] = np.nan
        subset = (adata.obs.Region == region) & (adata.obs.Weeks == weeks)
        adata.obs.loc[cells, "temp"] = adata.obs.loc[cells, "subclass"]
        adata.obs.loc[subset, "temp"] = adata.obs.loc[subset, "subclass"]
        adata.obs["temp"] = pd.Categorical(
            list(adata.obs["temp"]),
            categories=[
                "RGC",
                "Cone",
                "HC",
                "AC",
                "Rod",
                "BC",
            ],
        )
        sc.pl.umap(
            adata,
            color="temp",
            size=20,
            title="",
            frameon=False,
            legend_loc="None",
            palette={
                "Rod": "#17becf",
                "BC": "#bcbd22",
                "RGC": "#d62728",
                "Cone": "#e377c2",
                "HC": "#2ca02c",
                "AC": "#8c564b",
            },
        )
        fig = plt.gcf()
        fig.set_size_inches(4, 4)
        plt.savefig(
            "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_"
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
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Peripheral_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Peripheral_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Peripheral_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Peripheral_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Peripheral_PCW23.png",
    ],
    [
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Macula_PCW10.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Macula_PCW13.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Macula_PCW15.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Macula_PCW19.png",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_Macula_PCW23.png",
    ],
]

merged_image = merge_images(file_paths)
im = Image.fromarray(merged_image)
im.save(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure3/overall_umap_by_subclass_PCW.png",
    dpi=(300, 300),
)
