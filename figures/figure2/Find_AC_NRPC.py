#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sc.set_figure_params(scanpy=True, dpi=200)
colors = [
    "#440154FF",
    "#450558FF",
    "#46085CFF",
    "#470D60FF",
    "#471063FF",
    "#481467FF",
    "#481769FF",
    "#481B6DFF",
    "#481E70FF",
    "#482173FF",
    "#482576FF",
    "#482878FF",
    "#472C7AFF",
    "#472F7CFF",
    "#46327EFF",
    "#453581FF",
    "#453882FF",
    "#443B84FF",
    "#433E85FF",
    "#424186FF",
    "#404587FF",
    "#3F4788FF",
    "#3E4A89FF",
    "#3D4D8AFF",
    "#3C508BFF",
    "#3B528BFF",
    "#39558CFF",
    "#38598CFF",
    "#375B8DFF",
    "#355E8DFF",
    "#34608DFF",
    "#33638DFF",
    "#32658EFF",
    "#31688EFF",
    "#2F6B8EFF",
    "#2E6D8EFF",
    "#2D708EFF",
    "#2C718EFF",
    "#2B748EFF",
    "#2A768EFF",
    "#29798EFF",
    "#287C8EFF",
    "#277E8EFF",
    "#26818EFF",
    "#26828EFF",
    "#25858EFF",
    "#24878EFF",
    "#238A8DFF",
    "#228D8DFF",
    "#218F8DFF",
    "#20928CFF",
    "#20938CFF",
    "#1F968BFF",
    "#1F998AFF",
    "#1E9B8AFF",
    "#1F9E89FF",
    "#1FA088FF",
    "#1FA287FF",
    "#20A486FF",
    "#22A785FF",
    "#24AA83FF",
    "#25AC82FF",
    "#28AE80FF",
    "#2BB07FFF",
    "#2EB37CFF",
    "#31B67BFF",
    "#35B779FF",
    "#39BA76FF",
    "#3DBC74FF",
    "#41BE71FF",
    "#47C06FFF",
    "#4CC26CFF",
    "#51C56AFF",
    "#56C667FF",
    "#5BC863FF",
    "#61CA60FF",
    "#67CC5CFF",
    "#6DCD59FF",
    "#73D056FF",
    "#78D152FF",
    "#7FD34EFF",
    "#85D54AFF",
    "#8CD646FF",
    "#92D741FF",
    "#99D83DFF",
    "#A0DA39FF",
    "#A7DB35FF",
    "#ADDC30FF",
    "#B4DE2CFF",
    "#BBDE28FF",
    "#C2DF23FF",
    "#C9E020FF",
    "#D0E11CFF",
    "#D7E219FF",
    "#DDE318FF",
    "#E4E419FF",
    "#EBE51AFF",
    "#F1E51DFF",
    "#F7E620FF",
    # "#FDE725FF",
    "#5A5A5A",
]


# In[ ]:


adata = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/merged_h5ad/merged_raw_filtered_umap_10000_woadult_MG.h5ad"
)
NRPC = scv.read(
    "/storage/singlecell/zz4/fetal_snakemake/results/NRPC_fate/NRPC_fate.h5ad"
)
adata.obs["scpred_prediction"] = adata.obs.majorclass.replace(
    {
        "NRPC": "RPC",
        "PRPC": "RPC",
    }
)


# In[ ]:


sc.pp.normalize_total(adata, exclude_highly_expressed=True)
sc.pp.log1p(adata)


# In[ ]:


sc.pl.umap(NRPC, color="leiden", legend_loc="on data")


# # AC

# In[ ]:


sc.pl.umap(NRPC, color="AC")


# In[ ]:


for clusters in set(NRPC.obs.leiden):
    print(clusters)
    clusters = [clusters]
    adata.obs["temp"] = np.nan
    adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
        NRPC.obs.leiden.isin(clusters)
    ].AC
    adata.obs.loc[adata.obs.scpred_prediction == "AC", "temp"] = 1.01
    adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
    sc.pl.umap(
        adata,
        color="temp",
        palette=colors,
        legend_loc=None,
        frameon=False,
        title="",
    )


# In[ ]:


clusters = ["2", "13"]
print(clusters)
adata.obs["temp"] = np.nan
adata.obs.loc[NRPC.obs[NRPC.obs.leiden.isin(clusters)].index, "temp"] = NRPC.obs[
    NRPC.obs.leiden.isin(clusters)
].AC
adata.obs.loc[adata.obs.scpred_prediction == "AC", "temp"] = 1.01
adata.obs["temp"] = pd.cut(adata.obs["temp"], 100)
sc.pl.umap(
    adata, color="temp", palette=colors, legend_loc=None, frameon=False, title="",size = 20
)
fig = plt.gcf()
fig.set_size_inches(10, 10)
plt.savefig(
    "/storage/singlecell/zz4/fetal_snakemake/figures/figure2/find_AC_NRPC.svg",
    dpi=600,
    bbox_inches="tight",
)