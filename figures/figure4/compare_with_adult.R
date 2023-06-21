library(ggplot2)
output_dir = "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/"
coef <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/All_fit_coefs_region_models.csv")
coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral"),
    ]
genes <- c("YBX1", "PCP4", "PCAT4", "RP11-39E3.3", "TOMM7", "NACA", "FAU",
    "PFDN5", "EEF1A1", "TPT1", "UQCRB")
log2_fc <- c(2.372385326, 2.152899328, 2.037409037, 1.824192639, 1.735320302,
    1.712663892, 1.700248998, 1.683810323, 1.648282857, 1.566198246, 1.557601784)
df <- data.frame(genes, log2_fc)
rownames(df) <- df$genes
coef <- coef[coef$gene_short_name %in% genes, ]
coef$log2_fc <- df[coef$gene_short_name, "log2_fc"]
res <- coef

coef <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/All_fit_coefs_region_models.csv")
coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral"),
    ]
genes <- c("PLCD4", "DLG2", "MPP4", "LINC00969", "RAB41", "HSPB1", "CRYAB",
    "COL4A3", "MUM1", "MLXIP", "BCO2", "TF")
log2_fc <- c(-1.696328689, -1.705864553, -1.711252275, -1.720516862, -1.763088754,
    -1.879561396, -1.886103624, -1.899945676, -1.916313685, -1.999715663,
    -2.283023917, -2.809639996)
df <- data.frame(genes, log2_fc)
rownames(df) <- df$genes
coef <- coef[coef$gene_short_name %in% genes, ]
coef$log2_fc <- df[coef$gene_short_name, "log2_fc"]

res <- rbind(res, coef)
res <- res[order(res$log2_fc, decreasing = T), ]
res$gene_short_name <- factor(res$gene_short_name, levels = res$gene_short_name)
res$binary = "Macula enriched genes in adult samples"
res$binary[11:20] = "Peripheral enriched genes in adult samples"

svg(filename = paste(output_dir, "normalized_effect_barplot.svg", sep = ""),
    width = 15, height = 10)
ggplot(res, aes(x = res$gene_short_name, y = -res$normalized_effect, fill = binary)) +
    geom_bar(stat = "identity") + coord_flip() + ylab("LogFC") + xlab("Gene ID") +
    theme(legend.position = "bottom") + guides(fill = guide_legend(title = "")) +
    theme(text = element_text(size = 25))
dev.off()


svg(filename = paste(output_dir, "log2_fc_barplot.svg", sep = ""), width = 15,
    height = 10)
ggplot(res, aes(x = res$gene_short_name, y = res$log2_fc, fill = binary)) +
    geom_bar(stat = "identity") + coord_flip() + ylab("LogFC") + xlab("Gene ID") +
    theme(legend.position = "bottom") + guides(fill = guide_legend(title = "")) +
    theme(text = element_text(size = 25))
dev.off()