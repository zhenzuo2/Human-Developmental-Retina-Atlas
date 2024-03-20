import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import random
import magic
import seaborn as sns
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans


def even_density_downsampling_fast(coords, num_downsampled):
    """
    Even density downsampling of 2D coordinates using k-means clustering.

    Parameters:
    - coords: List of 2D coordinates [(x1, y1), (x2, y2), ...]
    - num_downsampled: Number of downsampled points to be selected.

    Returns:
    - List of selected indices representing the downsampled points.
    """

    # Convert coordinates to NumPy array
    data = np.array(coords)

    # Use k-means clustering to partition the data
    kmeans = KMeans(n_clusters=num_downsampled, random_state=42)
    kmeans.fit(data)

    # Get the cluster centers
    cluster_centers = kmeans.cluster_centers_

    # Find the closest point to each cluster center
    selected_indices = []
    for center in cluster_centers:
        closest_index = np.argmin(np.linalg.norm(data - center, axis=1))
        selected_indices.append(closest_index)

    return selected_indices


def remove_outlier_dots_umap(coordinates, threshold=2):
    # Compute the distances to the k nearest neighbors
    k = min(5, len(coordinates) - 1)  # Ensure k is smaller than the number of samples
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(coordinates)
    distances, _ = nbrs.kneighbors(coordinates)
    # Compute the median distance for each dot
    median_distances = np.median(distances[:, 1:], axis=1)
    # Identify outliers based on the median distance
    outlier_mask = median_distances > threshold
    # Return boolean array indicating if each dot was kept or removed
    return ~outlier_mask


def smooth(adata, layer):
    magic_op = magic.MAGIC()
    magic_op.set_params(n_jobs=10)
    emt_magic = magic_op.fit_transform(adata.layers[layer], genes="all_genes")
    emt_magic = magic_op.transform(genes="all_genes")
    adata.layers[layer + "_s"] = emt_magic
    return adata


def plot_hull(adata, adata_, layer, color):
    plt.clf()
    fig = plt.gcf()
    sns.kdeplot(
        x=[x[0] for x in adata_.obsm["X_umap"]],
        y=[x[1] for x in adata_.obsm["X_umap"]],
        weights=[x[0] for x in adata_[:, "CLU"].layers[layer].toarray()],
        cmap=color,
        shade=True,
    )
    # Given dot coordinates (x, y)
    dots = adata.obsm["X_umap"]
    # Calculate the convex hull
    hull = ConvexHull(dots)
    # Extract the vertices of the convex hull
    vertices = [dots[x] for x in hull.vertices.astype(int)]
    # Append the first vertex again to complete the polygon
    vertices = np.append(vertices, [vertices[0]], axis=0)
    # Extract the x and y coordinates separately
    x = vertices[:, 0]
    y = vertices[:, 1]
    # Plot the convex polygon
    plt.plot(x, y, "k-.")
    # Add labels and title
    plt.xlabel("")
    plt.ylabel("")
    plt.title("")
    plt.xticks([])
    plt.yticks([])

    # Display the plot
    fig.set_size_inches(6, 6)
    sns.despine(bottom=True, left=True)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/CLU_"
        + layer
        + "_ConvexHull.tiff",
        dpi=600,
        transparent=True,
        bbox_inches="tight",
    )


adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/multivelo_recover_dynamics_results/RPC_MG.h5ad"
)

adata = adata[adata.obs.majorclass.isin(["PRPC", "MG"])]

outlier_mask = remove_outlier_dots_umap(adata.obsm["X_umap"], 0.1)
adata = adata[outlier_mask]

adata = smooth(adata, "ATAC")
adata = smooth(adata, "unspliced")
adata = smooth(adata, "spliced")

sc.pl.umap(adata, color="CLU", layer="ATAC_s", vmax=2.5, frameon=False)
fig = plt.gcf()
plt.title("")
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/CLU_ATAC_s.tiff",
    dpi=600,
)

sc.pl.umap(adata, color="CLU", layer="unspliced_s", vmax=0.15, frameon=False)
fig = plt.gcf()
plt.title("")
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/CLU_unspliced_s.tiff",
    dpi=600,
)

sc.pl.umap(adata, color="CLU", layer="spliced_s", vmax=5, frameon=False)
fig = plt.gcf()
plt.title("")
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/chentemp/zz4/adult_dev_compare/figures/Figure5/CLU_spliced_s.tiff",
    dpi=600,
)

coordinates = adata.obsm["X_umap"]
num_downsampled = 5000
selected_indices = even_density_downsampling_fast(coordinates, num_downsampled)
adata_ = adata[adata.obs.index[selected_indices]]

plot_hull(adata, adata_, "ATAC", "Greys")
plot_hull(adata, adata_, "unspliced", "Greys")
plot_hull(adata, adata_, "spliced", "Greys")
