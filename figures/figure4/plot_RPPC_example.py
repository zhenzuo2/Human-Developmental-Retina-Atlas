import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import sys
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.spatial import ConvexHull
from sklearn.neighbors import NearestNeighbors
matplotlib.rcParams.update({'font.size': 22})

adata = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC_MG_NRPC.h5ad")
adata = adata[adata.obs.majorclass=="PRPC"]

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

outlier_mask=remove_outlier_dots_umap(adata.obsm["X_umap"],0.1)
adata = adata[outlier_mask]

mv.velocity_graph(adata)
mv.latent_time(adata)

scv.pl.scatter(adata, color='NFIA', size=20,layer="ATAC",vmax = 1.7)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_ATAC.svg",
    dpi=600,
)

scv.pl.scatter(adata, color='NFIA', size=20,layer="Ms",vmax = 0.3)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_Ms.svg",
    dpi=600,
)

scv.pl.scatter(adata, color='NFIA', size=20,layer="Mu",vmax = 7)
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_Mu.svg",
    dpi=600,
)

gene_list = ["NFIA"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(9, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_state.svg",
    dpi=600,
)

gene_list = ["CBX1"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(9, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/CBX1_state.svg",
    dpi=600,
)

gene_list = ["RPL10"]
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
mv.dynamic_plot(adata, gene_list, color_by='state', axis_on=False, frame_on=False)
fig = plt.gcf()
fig.set_size_inches(9, 3)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/RPL10_state.svg",
    dpi=600,
)
########################################################################################################


# Create a joint kernel density plot using Seaborn
sns.kdeplot(x = [x[0] for x in adata.obsm["X_umap"]], y= [x[1] for x in adata.obsm["X_umap"]], weights = [x[0] for x in adata[:,"NFIA"].layers["ATAC"].toarray()],cmap="Blues", shade=True)
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
plt.plot(x, y, 'k-.')

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("NFIA ATAC")

# Display the plot
fig = plt.gcf()
fig.set_size_inches(5, 5)
sns.despine(bottom = True, left = True)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_ATAC_ConvexHull.svg",
    dpi=600,
)

# Create a joint kernel density plot using Seaborn
sns.kdeplot(x = [x[0] for x in adata.obsm["X_umap"]], y= [x[1] for x in adata.obsm["X_umap"]], weights = [x[0] for x in adata[:,"NFIA"].layers["Ms"].toarray()],cmap="Greens", shade=True)
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
plt.plot(x, y, 'k-.')

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("NFIA spliced")

# Display the plot
fig = plt.gcf()
fig.set_size_inches(5, 5)
sns.despine(bottom = True, left = True)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_spliced_ConvexHull.svg",
    dpi=600,
)

# Create a joint kernel density plot using Seaborn
sns.kdeplot(x = [x[0] for x in adata.obsm["X_umap"]], y= [x[1] for x in adata.obsm["X_umap"]], weights = [x[0] for x in adata[:,"NFIA"].layers["Mu"].toarray()],cmap="Reds", shade=True)
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
plt.plot(x, y, 'k-.')

# Add labels and title
plt.xlabel("")
plt.ylabel("")
plt.title("NFIA unspliced")

# Display the plot
fig = plt.gcf()
sns.despine(bottom = True, left = True)
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/NFIA_unspliced_ConvexHull.svg",
    dpi=600,
)

###
mv.pie_summary(adata)
# Display the plot
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/PRPCpie_summary.svg",
    dpi=600,
    bbox_inches='tight'
)

################################################################################

PRPC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/PRPC.h5ad")
NRPC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/NRPC.h5ad")
AC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/AC.h5ad")
BC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/BC.h5ad")
Rod = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/Rod.h5ad")
Cone = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/Cone.h5ad")
RGC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/RGC.h5ad")
HC = sc.read("/storage/singlecell/zz4/fetal_snakemake/results/multivelo_recover_dynamics_results/HC.h5ad")

def get_summary(adata):
    t_sw = adata[:, adata.var['velo_s_genes']].var[['fit_t_sw1', 'fit_t_sw2', 'fit_t_sw3']].copy()
    t_sw = t_sw.mask(t_sw > 20, 20)
    t_sw = t_sw.mask(t_sw < 0)
    t_sw['interval 1'] = t_sw['fit_t_sw1']
    t_sw['t_sw2 - t_sw1'] = t_sw['fit_t_sw2'] - t_sw['fit_t_sw1']
    t_sw['t_sw3 - t_sw2'] = t_sw['fit_t_sw3'] - t_sw['fit_t_sw2']
    t_sw['20 - t_sw3'] = 20 - t_sw['fit_t_sw3']
    t_sw = t_sw.mask(t_sw <= 0)
    t_sw = t_sw.mask(t_sw > 20)
    t_sw.columns = pd.Index(['time 1', 'time 2', 'time 3', 'primed',
                             'coupled-on', 'decoupled', 'coupled-off'])
    t_sw = t_sw[['primed', 'coupled-on', 'decoupled', 'coupled-off']]
    t_sw = t_sw / 20
    a = t_sw['primed'][~np.isnan(t_sw['primed'])]
    a = np.median(a)
    b = t_sw['coupled-on'][~np.isnan(t_sw['coupled-on'])]
    b = np.median(b)
    c = t_sw['decoupled'][~np.isnan(t_sw['decoupled'])]
    c = np.median(c)
    d = t_sw['coupled-off'][~np.isnan(t_sw['coupled-off'])]
    d = np.median(d)
    return([a,b,c,d])


# Generate random data for three time points
time_points = ['PRPC', 'NRPC',"RGC","Cone","HC",'AC',"Rod", "BC",]
observations = [get_summary(PRPC),
                get_summary(NRPC),
                get_summary(RGC),
                get_summary(Cone),
                get_summary(HC),
                get_summary(AC),
                get_summary(Rod),
                get_summary(BC),
               ]
observations = observations / np.sum(observations, axis=1)[:, np.newaxis]

# Set up colors for the stacked bars
colors = ['#DC3912', '##FF9900', '#109618', '##3366CC']

# Plot the stacked bars

bottom = np.zeros(len(time_points))
plt.clf()

for i, obs in enumerate(observations.T):
    plt.bar(time_points, obs, bottom=bottom, color=colors[i])
    bottom += obs

# Add labels and title
plt.xlabel('')
plt.ylabel('')
plt.title('')
# Show the legend
plt.legend(['Primed', 'Coupled-on', 'Decoupled', 'Coupled-off'],bbox_to_anchor=(1.04, 1), loc="upper left")
fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks(
    ticks=[0, 1, 2, 3,4,5,6,7],
    labels=['PRPC', 'NRPC',"RGC","Cone","HC",'AC',"Rod", "BC",],
    rotation=45,
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/switch_time_summary.svg",
    dpi=600,
    bbox_inches='tight'
)

def get_summary(adata):
    genes = adata.var_names
    fit_model = adata[:, (adata.var['fit_direction'] == 'complete') &
                      np.isin(adata.var_names, genes)].var['fit_model'].values
    fit_direction = adata[:, genes].var['fit_direction'].values
    data = [np.sum(fit_direction == 'on'), np.sum(fit_direction == 'off'),
            np.sum(fit_model == 1), np.sum(fit_model == 2)]
    index = ['induction', 'repression', 'Model 1', 'Model 2']
    index = [x for i, x in enumerate(index) if data[i] > 0]
    data = [x for x in data if x > 0]
    df = pd.DataFrame({'data': data}, index=index)
    return df

observations = [get_summary(PRPC),
                get_summary(NRPC),
                get_summary(RGC),
                get_summary(Cone),
                get_summary(HC),
                get_summary(AC),
                get_summary(Rod),
                get_summary(BC),
               ]
res = np.zeros((8, 4))
for i in range(8):
    j=0
    for x in ['induction', 'repression', 'Model 1', 'Model 2']:
        if x in observations[i].index:
            res[i,j] =  observations[i].loc[x,"data"]
        j = j+1


res = res / np.sum(res, axis=1)[:, np.newaxis]

# Set up colors for the stacked bars
colors = ['#0066CC', '#99CC00', '#FF6699', '#FFCC00']

# Plot the stacked bars

bottom = np.zeros(len(time_points))
plt.clf()
for i, obs in enumerate(res.T):
    plt.bar(time_points, obs, bottom=bottom, color=colors[i])
    bottom += obs

# Add labels and title
plt.xlabel('')
plt.ylabel('')
plt.title('')

# Show the legend
plt.legend(['Induction', 'Repression', 'Model 1', 'Model 2'],labelcolor = colors,bbox_to_anchor=(1.04, 1), loc="upper left")

fig = plt.gcf()
fig.set_size_inches(5, 5)
plt.xticks(
    ticks=[0, 1, 2, 3,4,5,6,7],
    labels=['PRPC', 'NRPC',"RGC","Cone","HC",'AC',"Rod", "BC",],
    rotation=45,
)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/model_summary.svg",
    dpi=600,
    bbox_inches='tight'
)
