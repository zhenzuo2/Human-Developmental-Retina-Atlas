library(ggplot2)
output_dir = "/storage/singlecell/zz4/fetal_snakemake/figures/figure4/"
coef <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/All_fit_coefs_region_days_models.csv")
coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral"),
]
genes <- c('PCAT4', 'NACA', 'FAU', 'EEF1A1', 'TPT1', 'PDE6H', 'MT-CO2', 'MT-CO3',
       'SNAP25', 'MAP1B')
log2_fc <- c(2.03740904, 1.71266389, 1.700249  , 1.64828286, 1.56619825,
       1.40848186, 1.06995073, 1.05443919, 1.04698432, 1.03871918)
df <- data.frame(genes, log2_fc)
rownames(df) <- df$genes
coef <- coef[coef$gene_short_name %in% genes, ]
coef$log2_fc <- df[coef$gene_short_name, "log2_fc"]
coef$normalized_effect <-c(
0.16434465, -0.17517452, -0.3409974 , -0.4658646 , -0.3368352 ,
        2.4855852 , -0.20533793, -0.50852257, -0.2591621 ,  0.37719795)
res <- coef

coef <- read.csv("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/All_fit_coefs_region_days_models.csv")
coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral"),
]
genes <- c('DST', 'ARHGAP21', 'RBM6', 'HIBCH', 'ODF2L', 'NOS1AP', 'SRGAP3',
       'SORBS2', 'DLG2', 'MPP4')
log2_fc <- c(-1.39009953, -1.39368146, -1.4146402 , -1.47889365, -1.53206764,
       -1.54879592, -1.65002095, -1.69446377, -1.70586455, -1.71125228)
df <- data.frame(genes, log2_fc)
rownames(df) <- df$genes
coef <- coef[coef$gene_short_name %in% genes, ]
coef$log2_fc <- df[coef$gene_short_name, "log2_fc"]
coef$normalized_effect <-c(
  	0.2472441 , 0.0691049 , 0.35834637, 0.41975757, 0.7226046 ,
       2.2553082 , 0.20061453, 0.67848617, 1.6413188 , 1.338123
)
res <- rbind(res, coef)
res <- res[order(res$log2_fc, decreasing = T), ]
res$gene_short_name <- factor(res$gene_short_name, levels = res$gene_short_name)
res$binary = "Macula enriched genes in adult samples"
res$binary[11:20] = "Peripheral enriched genes in adult samples"

svg(filename = paste(output_dir, "normalized_effect_barplot.svg", sep = ""),
    width = 15, height = 10)
ggplot(res, aes(x = res$gene_short_name, y = res$normalized_effect, fill = binary)) +
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