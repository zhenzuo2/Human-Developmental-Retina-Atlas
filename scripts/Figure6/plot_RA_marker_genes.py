import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
markers = [
    "NFKB1",
    "FABP5",
    "VAX2",
    "TBX5",
    "EFNB1",
    "EFNB2",
    "EPHB2",
    "EPHB3",
    "DHRS3",
    "RBP4",
    "RDH10",
    "FGF8",
    "FGF9",
    "STRA6",
    "CYP26A1",
    "CYP26B1",
    "CYP26C1",
    "RXRG",
    "RXRA",
    "RXRB",
    "RARA",
    "RARB",
    "RARG",
    "ALDH1A1",
    "ALDH1A2",
    "ALDH1A3",
]

for gene in markers:
    vmax = np.quantile(adata[:, gene].X.toarray(), q=0.999)
    plt.clf()
    sc.pl.umap(adata[adata.obs.Region == "Macula"], color=gene, vmax=vmax, size=30)
    plt.xlabel("")
    plt.ylabel("")
    plt.title(gene + " (Macula)", fontsize=25)
    plt.legend(loc='lower left')
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_Macula.tiff",
        dpi=600,
        bbox_inches="tight",
        transparent=True
    )

    sc.pl.umap(adata[adata.obs.Region == "Peripheral"], color=gene, vmax=vmax, size=30)
    plt.xlabel("")
    plt.ylabel("")
    plt.title(gene+ " (Peripheral)", fontsize=25)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_Peripheral.tiff",
        dpi=600,
        bbox_inches="tight",
        transparent=True
    )

    plt.clf()
    #adata = adata[adata.obs.Region !="Whole Eye"]
    df = sc.get.obs_df(adata,  [gene, "Region"])
    df = df.set_index('Region').stack().reset_index()
    df.columns = ['Region', 'gene', 'value']
    df['NonZero'] = df['value'] != 0
    # Group by 'Region' and 'Gene', then calculate the percentage of non-zero values
    result = df.groupby(['Region', 'gene'])['NonZero'].mean().reset_index()
    df = pd.DataFrame(result)
    # Create subplots for each gene
    genes = df['gene'].unique()
    
    # Define colors for each region
    region_colors = {"Whole Eye": "#ffdd05", "Macula" : "#F8766D", "Peripheral": "#00BFC4"}
    # Iterate over genes and create bar plots
    gene_df = df[df['gene'] == gene]
    for region in gene_df['Region'].unique():
        region_df = gene_df[gene_df['Region'] == region]
        plt.bar(region_df['Region'], region_df['NonZero'] * 100, color=region_colors.get(region), label=region)
    plt.title(gene, fontsize=25)
    plt.ylabel('Cells (Percent)')
    plt.xlabel('')
    plt.legend()
    plt.xticks(fontsize=15)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_barplot.tiff",
        dpi=600,
        bbox_inches="tight",
        transparent=True
    )


def merge_tiffs(images, output_path, padding_color=(0, 0, 0, 0)):
    # Open all TIFF images and get their sizes
    image_objects = [Image.open(image).convert("RGBA") for image in images]
    image_sizes = [image.size for image in image_objects]

    # Find the maximum height and total width
    max_height = max([height for _, height in image_sizes])
    total_width = sum([width for width, _ in image_sizes])

    # Create a new image with the maximum height and total width
    output_image = Image.new("RGBA", (total_width, max_height), padding_color)

    # Paste each image into the new image with appropriate offset
    current_width = 0
    for image, size in zip(image_objects, image_sizes):
        offset = (current_width, (max_height - size[1]) // 2)
        output_image.paste(image, offset)
        current_width += size[0]

    # Save the final image
    output_image.save(output_path)

for gene in markers:
    input_paths = [
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_Macula.tiff",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_Peripheral.tiff",
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/RA_genes/"
        + gene
        + "_barplot.tiff",

    ]
    output_path = (
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/" + gene + ".tiff"
    )
    merge_tiffs(input_paths, output_path)
