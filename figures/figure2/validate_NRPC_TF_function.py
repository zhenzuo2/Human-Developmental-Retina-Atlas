import pandas as pd
import scanpy as sc
import numpy as np

adata = sc.read("merged_raw_filtered_umap_10000_woadult_MG.h5ad")
adata = adata[adata.obs.majorclass.isin(["AC", "BC", "Cone", "HC", "RGC", "Rod"])]
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
meta = pd.read_csv("NRPC_grn_plot.csv")
df = pd.read_csv("Annotated_NRPC_modules_meta.csv")
df["tf_ge"] = ""
df["tf_score"] = ""
for gene in set(meta.name):
    print(gene)
    adata.obs[gene] = adata[:, gene].X.toarray()
    temp = adata.obs.groupby("majorclass").mean(numeric_only=True)[gene]
    if np.sum(temp / np.sum(temp) > 0.33) > 0:
        df.loc[df.tf == gene, "tf_ge"] = "_".join(
            temp[temp / np.sum(temp) > 0.33].index.values
        )
    else:
        df.loc[df.tf == gene, "tf_ge"] = temp.idxmax()
    tf_score = pd.DataFrame(index=list(set(adata.obs.majorclass)))
    tf_score["values"] = 0
    for celltype in set(adata.obs.majorclass):
        lista = [
            np.mean(
                adata[
                    adata.obs.majorclass == celltype, df.loc[df.tf == gene, "target"]
                ].X,
                axis=0,
            )
        ]
        listb = df.loc[df.tf == gene, "estimate"]
        tf_score.loc[celltype, "values"] = np.sum([a * b for a, b in zip(lista, listb)])
    if np.sum(tf_score.values / np.sum(tf_score.values) > 0.33) > 0:
        df.loc[df.tf == gene, "tf_score"] = "_".join(
            tf_score.loc[(tf_score / np.sum(tf_score) > 0.33).values, :].index.values
        )
    else:
        df.loc[df.tf == gene, "tf_score"] = tf_score.idxmax()[0]

res = pd.DataFrame.drop_duplicates(df[["tf", "tf_ge", "tf_score"]]).reset_index(
    drop=True
)
res = res.loc[res.tf.isin(meta.name), :]
res = res.reset_index(drop=True)
pd.set_option("display.max_rows", None)
res.index = res.tf.values
val = pd.read_csv(
    "/Users/zhenzuo/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/TOP TF during fetal retina development/table/table.csv"
)
val.index = val.Name.values
val = val.loc[[x for x in val.index if x in res.index],]
val = pd.concat([res, val.Spec_pos], axis=1)
val["Spec_pos"] = val["Spec_pos"].astype(str)
na = []
same_ge = []
same_as = []
diff = []
for i in range(val.shape[0]):
    set1 = val.tf_ge[i].split("_")
    set2 = val.tf_score[i].split("_")
    if val.Spec_pos[i] == "nan":
        na = na + [val.index[i]]
        continue
    else:
        set3 = val.Spec_pos[i].split(",")
    if len([x for x in set1 if x in set3]) > 0:
        same_ge = same_ge + [val.index[i]]
    if len([x for x in set2 if x in set3]) > 0:
        same_as = same_as + [val.index[i]]
    if len([x for x in set1 + set2 if x in set3]) == 0:
        diff = diff + [val.index[i]]
(" ").join(same_ge)
(" ").join(same_as)
(" ").join(na)
