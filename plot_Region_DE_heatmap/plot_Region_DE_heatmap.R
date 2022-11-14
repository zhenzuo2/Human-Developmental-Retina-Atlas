library(ggplot2)

output_dir = c("/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_",
               "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_",
               "/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/",
               "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/")

cm_files = c("All_compare_mod_region.csv", "Early_compare_mod_region.csv",
             "Late_compare_mod_region.csv")

coef_files = c("All_fit_coefs_region_models.csv", "Early_fit_coefs_region_models.csv",
               "Late_fit_coefs_region_models.csv")

mode = c("All", "Early", "Late")

for (out in output_dir) {
  for (i in 1:3) {
    coef <- read.csv(paste(out, coef_files[i], sep = ""))
    coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral") &
                   (coef$q_value < 0.01) & (abs(coef$normalized_effect) >= 0.5),
    ]
    cm <- read.csv(paste(out, cm_files[i], sep = ""))
    cm <- cm[cm$q_value < 0.01, ]
    write.table(x = data.frame(intersect(coef$gene_id, cm$gene_short_name)),
                file = paste(out, mode[i], "_Region_DE_filtered_gene_list.csv",
                             sep = ""), row.names = F, col.names = F, quote = F)
    print(paste(out, mode[i], "_Region_DE_filtered_gene_list.csv",
                sep = ""))
  }
}

files <- c("/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/PRPC_Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/MG_monocle3_DE_analysis/NRPC_Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/AC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/BC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Cone_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/Rod_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/HC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_bash/results/RGC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv")

names <- c("RPC_All", "RPC_Early", "RPC_Late", "PRPC_All", "PRPC_Early",
           "PRPC_Late", "NRPC_All", "NRPC_Early", "NRPC_Late", "AC_All", "AC_Early",
           "AC_Late", "BC_All", "BC_Early", "BC_Late", "Cone_All", "Cone_Early",
           "Cone_Late", "Rod_All", "Rod_Early", "Rod_Late", "HC_All", "HC_Early",
           "HC_Late", "RGC_All", "RGC_Early", "RGC_Late")

n <- length(files)
gene_list <- data.frame(Name = character(), Set = character())
for (i in 1:n) {
  temp <- read.csv(files[i])
  colnames(temp) <- "Name"
  temp$Set <- names[i]
  gene_list <- rbind(gene_list, temp)
}

mat <- matrix(data = 0, nrow = n, ncol = n)
rownames(mat) <- names
colnames(mat) <- names

for (i in 1:n) {
  for (j in 1:n) {
    mat[i, j] <- length(intersect(gene_list[gene_list$Set == names[i],
                                            "Name"], gene_list[gene_list$Set == names[j], "Name"]))
  }
}
melted_mat <- melt(mat)


svg("/storage/singlecell/zz4/fetal_bash/figures/Region_DE_Heatmap/samples.svg",width = 10, height = 10)
ggplot(data = melted_mat, aes(x = Var1, y = Var2, fill = value)) + geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradientn(colours = heat.colors(10)) + xlab("") + ylab("") +
  geom_text(aes(x = Var1, y = Var2, label = value), color = "black",
            size = 3)
dev.off()
