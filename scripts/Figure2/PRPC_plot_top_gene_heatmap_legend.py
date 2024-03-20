import matplotlib.pyplot as plt
import numpy as np
output_file_path = "/storage/chentemp/zz4/adult_dev_compare/figures/Figure2/"
# Define a colormap and normalization
cmap = plt.cm.RdBu_r
norm = plt.Normalize(-1, 1)

# Create a figure and axis
fig, ax = plt.subplots(figsize=(6, 1))

# Create a ScalarMappable to map the normalization to colors
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add a colorbar to the figure
cbar = plt.colorbar(sm, orientation='horizontal', ax=ax)
cbar.set_label('z-score')

# Hide the axis
ax.set_axis_off()

# Display the figure
plt.savefig(
    output_file_path + "PRPC_plot_top_gene_heatmap_legend.tiff",
    bbox_inches="tight",
    transparent=True,
    dpi = 300
)
