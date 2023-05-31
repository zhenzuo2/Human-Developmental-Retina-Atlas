library("EnhancedVolcano")
res <- read.csv("/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_subtype_DEGs.csv")
rownames(res) <- res$names
EnhancedVolcano(res, lab = rownames(res), x = "logfoldchanges", y = "pvals_adj",
    pCutoff = 0.1, FCcutoff = 1, pointSize = 1, labSize = 6, colAlpha = 1,
    legendPosition = "none", max.overlaps = 15, title = "", subtitle = "",
    ylab = bquote(~-Log[10] ~ italic((adjustedP))))
ggsave("/storage/singlecell/zz4/fetal_snakemake/figures/figure3/AC_subtype_DEGs.svg",
    bg = "transparent", width = 4.2, height = 3, units = "in")
