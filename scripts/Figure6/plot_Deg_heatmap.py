import pandas as pd
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
import seaborn as sns
import csv

def is_valid_string(s):
    return not s.startswith("MT") and re.match("^[a-zA-Z0-9_]*$", s) is not None

def count_shared_elements(dict1, dict2):
    keys = dict1.keys()  # Assuming keys are the same in both dictionaries
    print(keys)
    # Initialize a matrix to store the count of shared elements
    matrix = [[0 for _ in keys] for _ in keys]
    # Iterate through each pair of keys
    for i, key1 in enumerate(keys):
        for j, key2 in enumerate(keys):
            # Count the number of shared elements between the lists
            shared_elements = len(set(dict1[key1]) & set(dict2[key2]))
            # Save the count in the matrix
            matrix[i][j] = shared_elements
    matrix = pd.DataFrame(matrix)
    matrix.index = ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]
    matrix.columns = ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]
    return matrix

def plot_clustered_heatmap(
    data, annot=True, cmap="coolwarm", fig_name=None
):
    # Plot clustered heatmap
    plt.clf()
    sns.heatmap(
        data,
        cmap=cmap,
        annot=annot,
        fmt=".0f",
        vmax = 200
    )
    # Show the plot
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.xticks(rotation=90) 
    plt.savefig(
        fig_name,
        dpi=300,
        bbox_inches="tight",
        transparent=True,
    )

def save_dict_to_csv(dictionary, filename):
    # Extract keys and values from the dictionary
    keys = list(dictionary.keys())
    values = list(dictionary.values())

    # Find the maximum length of any list in the values
    max_list_length = max(len(lst) for lst in values)

    # Create a list of lists to store rows for CSV
    rows = [[] for _ in range(max_list_length + 1)]  # +1 for header row

    # Add header row with keys
    rows[0] = keys

    # Add values to rows
    for i in range(max_list_length):
        for j in range(len(keys)):
            if i < len(values[j]):
                rows[i + 1].append(values[j][i])
            else:
                rows[i + 1].append('')  # Add empty string if the list is shorter

    # Write to CSV file
    with open(filename, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerows(rows)



P = {}
for x in ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]:
    print(x)
    df = pd.read_csv(
        "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/" + x + "_monocle3_DE_analysis/compare_mod_region.csv"
    )
    df = df.loc[(df.q_value < 0.01) & (df.num_cells_expressed > 1000), :]
    gene_list_1 = list(df.gene_short_name)
    df = pd.read_csv(
        "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/"
        + x
        + "_monocle3_DE_analysis/fit_coefs_region_models.csv"
    )
    df = df.loc[
        (df.term == "RegionPeripheral")
        & (df.q_value < 0.01)
        & (df.normalized_effect > 1),
    ]
    gene_list_2 = list(df.gene_short_name)
    gene_list = list(set(gene_list_1).intersection(gene_list_2))
    filtered_strings = [s for s in gene_list if is_valid_string(s)]
    # Print the result
    print(len(filtered_strings))
    P[x] = filtered_strings

M = {}
for x in ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]:
    print(x)
    df = pd.read_csv(
        "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/" + x + "_monocle3_DE_analysis/compare_mod_region.csv"
    )
    df = df.loc[(df.q_value < 0.01) & (df.num_cells_expressed > 1000), :]
    gene_list_1 = list(df.gene_short_name)

    df = pd.read_csv(
        "/storage/chentemp/zz4/adult_dev_compare/results/monocle3_DE_analysis/"
        + x
        + "_monocle3_DE_analysis/fit_coefs_region_models.csv"
    )
    df = df.loc[
        (df.term == "RegionPeripheral")
        & (df.q_value < 0.01)
        & (df.normalized_effect < -1),
    ]
    gene_list_2 = list(df.gene_short_name)
    gene_list = list(set(gene_list_1).intersection(gene_list_2))
    filtered_strings = [s for s in gene_list if is_valid_string(s)]

    # Print the result
    M[x] = filtered_strings

save_dict_to_csv(P, '/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/DEGs_P.csv')
save_dict_to_csv(M, '/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/DEGs_M.csv')

# Get the matrix of shared elements
result_matrix_P = count_shared_elements(P, P)
result_matrix_M = count_shared_elements(M, M)

labels = ["AC", "BC", "Cone", "HC", "MG", "NRPC", "PRPC", "RGC", "Rod"]

print(result_matrix_P)
# Plot the heatmap with custom axis labels
plot_clustered_heatmap(
    result_matrix_P,
    fig_name="/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/plot_Deg_heatmap_P.tiff",
)
plot_clustered_heatmap(
    result_matrix_M,
    fig_name="/storage/chentemp/zz4/adult_dev_compare/figures/Figure6/plot_Deg_heatmap_M.tiff",
)
