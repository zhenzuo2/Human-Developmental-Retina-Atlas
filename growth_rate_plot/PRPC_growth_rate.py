
import matplotlib.pyplot as plt
import scvelo as scv
import scanpy as sc
import cellrank as cr

# import CellRank kernels and estimators
from cellrank.external.kernels import WOTKernel
from cellrank.tl.kernels import ConnectivityKernel
from cellrank.tl.estimators import GPCCA

# set verbosity levels
cr.settings.verbosity = 2
scv.settings.verbosity = 3
adata_result = scv.read(
    "/storage/singlecell/zz4/fetal_bash/results/merged_h5ad/merged_raw_filtered_umap_10000_major_sub_class.h5ad"
)
adata_result

Time = pd.read_csv(
    "/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_results/PRPC.obs.csv"
)
adata_result = adata_result[Time["Unnamed: 0"]]

adata_result.obs["latent_time"] = Time["latent_time"].values

# figure settings
scv.settings.set_figure_params(
    "scvelo", dpi_save=400, dpi=80, transparent=True, fontsize=20, color_map="viridis"
)
wk = WOTKernel(adata_result, time_key="Days")
adata_result.var_names_make_unique()
adata_result.obs_names_make_unique()
wk.compute_initial_growth_rates(organism="human", key_added="growth_rate_init")


adata = scv.read("/storage/singlecell/zz4/fetal_bash/results/multivelo_recover_dynamics_run_umap_PRPC/adata_umap.h5ad")
adata.obs['growth_rate_init'] = adata_result[adata.obs.index].obs['growth_rate_init'].values
scv.pl.scatter(
    adata, color="growth_rate_init", legend_loc="right", basis="umap", s=20,save = "growth_rate_init.svg",title = "",color_map = "inferno"
)