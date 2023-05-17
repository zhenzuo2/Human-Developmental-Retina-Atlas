library("EnhancedVolcano")
res<- read.csv("/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_DE.csv")
rownames(res) <- res$names
svg("/storage/singlecell/zz4/fetal_snakemake/figures/figure2/AC_NRPC_DE.svg")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'logfoldchanges',
                y = 'pvals_adj',
                pCutoff = 10e-2,
                FCcutoff = 1.0,
                pointSize = 1.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = 'none',
                max.overlaps =15,
                title = "",
                subtitle = "",ylab = bquote(~-Log[10] ~ italic((adjustedP))))
dev.off()