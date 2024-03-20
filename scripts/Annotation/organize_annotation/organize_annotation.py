import scanpy as sc
import pandas as pd

adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult.h5ad"
)

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/AC_subtype_processed.csv"
)
AC = adata[df["Unnamed: 0"]]
AC.obs["majorclass"] = df.majorclass.values
AC.obs["subclass"] = df.subclass.values
AC.obs["celltype"] = df.celltype.values

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/BC_subtype_processed.csv"
)
BC = adata[df["Unnamed: 0"]]
BC.obs["majorclass"] = df.majorclass.values
BC.obs["subclass"] = df.subclass.values
BC.obs["celltype"] = df.celltype.values

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/HC_subtype_processed.csv"
)
HC = adata[df["Unnamed: 0"]]
HC.obs["majorclass"] = df.majorclass.values
HC.obs["subclass"] = df.subclass.values
HC.obs["celltype"] = df.celltype.values

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/RGC_subtype_processed.csv"
)
RGC = adata[df["Unnamed: 0"]]
RGC.obs["majorclass"] = df.majorclass.values
RGC.obs["subclass"] = df.subclass.values
RGC.obs["celltype"] = df.celltype.values

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Rod_subtype_processed.csv"
)
Rod = adata[df["Unnamed: 0"]]
Rod.obs["majorclass"] = df.majorclass.values
Rod.obs["subclass"] = df.subclass.values
Rod.obs["celltype"] = df.celltype.values

df = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/cell_annotation_results/subclass_annotation_result_csv/Cone_subtype_processed.csv"
)
Cone = adata[df["Unnamed: 0"]]
Cone.obs["majorclass"] = df.majorclass.values
Cone.obs["subclass"] = df.subclass.values
Cone.obs["celltype"] = df.celltype.values

MG = adata[adata.obs.majorclass == "MG"]
MG.obs["majorclass"] = "NRPC"
MG.obs["subclass"] = "NRPC"
MG.obs["celltype"] = "NRPC"

MG.obs.loc[MG.obs.leiden == "12", "majorclass"] = "MG"
MG.obs.loc[MG.obs.leiden == "12", "subclass"] = "MG"
MG.obs.loc[MG.obs.leiden == "12", "celltype"] = "MG"

MG.obs.loc[MG.obs.leiden.isin(["24", "1", "2", "23", "39"]), "majorclass"] = "PRPC"
MG.obs.loc[MG.obs.leiden.isin(["24", "1", "2", "23", "39"]), "subclass"] = "PRPC"
MG.obs.loc[MG.obs.leiden.isin(["24", "1", "2", "23", "39"]), "celltype"] = "PRPC"

temp = anndata.concat([AC, BC, RGC, Rod, Cone, HC, MG])

temp.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated_not_filter_ATAC.h5ad"
)

atac_cells = pd.read_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/ATAC_filtered_cells/ATAC_Filtered_cell.csv",
    sep=" ",
)
atac_filtered_cells = [string.replace("#", "_") for string in atac_cells.x]

# Find common cells after ATAC-seq and RNA-seq filtering
common_cells = set(temp.obs.index).intersection(atac_filtered_cells)
temp = temp[list(common_cells)]

temp.write(
    "/storage/chentemp/zz4/adult_dev_compare/results/merged_h5ad/merged_raw_filtered_umap_10000_wadult_annotated.h5ad"
)

for cell_type in list(set(temp.obs.majorclass)):
    temp[temp.obs.majorclass == cell_type].obs.to_csv(
        "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/"
        + cell_type
        + ".csv"
    )

temp[temp.obs.majorclass.isin(["PRPC", "NRPC", "MG"])].obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/RPC_MG.csv"
)
temp[temp.obs.majorclass.isin(["PRPC", "MG"])].obs.to_csv(
    "/storage/chentemp/zz4/adult_dev_compare/results/Annotation/PRPC_MG.csv"
)
