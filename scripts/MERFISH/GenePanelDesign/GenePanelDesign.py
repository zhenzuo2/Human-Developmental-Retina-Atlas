import scanpy as sc
import scmer
import csv
from scmer import UmapL1

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)

sc.pp.filter_genes(adata, min_cells=50)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.01, min_disp=0.2)
adata = adata[:, adata.var.highly_variable]

print(adata)

sc.pp.subsample(adata, n_obs=10000)

model_20 = UmapL1.tune(
    target_n_features=550,
    X=adata.X.toarray(),
    n_threads=10,
    use_gpu=True,
)
print(*adata.var_names[model_20.get_mask()])

my_list = adata.var_names[model_20.get_mask()]
csv_file = "/storage/chentemp/zz4/adult_dev_compare/results/Gene_panel/500.csv"

# Open the CSV file in write mode
with open(csv_file, "w", newline="") as file:
    # Create a CSV writer object
    writer = csv.writer(file)

    # Write each item in the list as a separate row
    for item in my_list:
        writer.writerow([item])

print(f"The list has been successfully saved to {csv_file}.")
