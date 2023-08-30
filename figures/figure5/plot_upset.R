library(ggplot2)
library(UpSetR)
for (region in c("Macula", "Peripheral")) {
  output_dir = c("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_",
                 #"/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/MG_",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/AC_monocle3_DE_analysis/",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/BC_monocle3_DE_analysis/",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Rod_monocle3_DE_analysis/",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/HC_monocle3_DE_analysis/",
                 "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/RGC_monocle3_DE_analysis/")
  
  cm_files = c("All_compare_mod_region.csv", "Early_compare_mod_region.csv",
               "Late_compare_mod_region.csv")
  
  coef_files = c("All_fit_coefs_region_models.csv", "Early_fit_coefs_region_models.csv",
                 "Late_fit_coefs_region_models.csv")
  
  mode = c("All", "Early", "Late")
  
  for (out in output_dir) {
    for (i in 1:3) {
      coef <- read.csv(paste(out, coef_files[i], sep = ""))
      if (region == "Macula") {
        coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral") &
                       (coef$q_value < 0.01/27) & (coef$normalized_effect <= -0.5),
        ]
      }
      if (region == "Peripheral") {
        coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral") &
                       (coef$q_value < 0.01/27) & (coef$normalized_effect >= 0.5),
        ]
      }
      cm <- read.csv(paste(out, cm_files[i], sep = ""))
      cm <- cm[cm$q_value < 0.01/27, ]
      write.table(x = data.frame(intersect(coef$gene_id, cm$gene_short_name)),
                  file = paste(out, mode[i], "_Region_DE_filtered_gene_list.csv",
                               sep = ""), row.names = F, col.names = F, quote = F)
      print(paste(out, mode[i], "_Region_DE_filtered_gene_list.csv",
                  sep = ""))
    }
  }
  
  files <- c("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_Late_Region_DE_filtered_gene_list.csv",
             #"/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_All_Region_DE_filtered_gene_list.csv",
             #"/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_Early_Region_DE_filtered_gene_list.csv",
             #"/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/MG_All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/MG_Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/MG_Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/AC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/AC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/AC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/BC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/BC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/BC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Cone_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Rod_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Rod_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/Rod_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/HC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/HC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/HC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/RGC_monocle3_DE_analysis/All_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/RGC_monocle3_DE_analysis/Early_Region_DE_filtered_gene_list.csv",
             "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/RGC_monocle3_DE_analysis/Late_Region_DE_filtered_gene_list.csv")
  
  names <- c("PRPC_All", "PRPC_Early", "PRPC_Late", "NRPC_All", "NRPC_Early",
             "NRPC_Late", #"RPC_All", "RPC_Early", "RPC_Late", 
             "MG_All", "MG_Early",
             "MG_Late", "AC_All", "AC_Early", "AC_Late", "BC_All", "BC_Early",
             "BC_Late", "Cone_All", "Cone_Early", "Cone_Late", "Rod_All", "Rod_Early",
             "Rod_Late", "HC_All", "HC_Early", "HC_Late", "RGC_All", "RGC_Early",
             "RGC_Late")
  
  n <- length(files)
  gene_list <- data.frame(Name = character(), Set = character())
  for (i in 1:n) {
    temp <- read.csv(files[i])
    colnames(temp) <- "Name"
    temp$Set <- names[i]
    gene_list <- rbind(gene_list, temp)
  }
  gene_list <- gene_list[!grepl("[^A-Za-z0-9]", gene_list$Name), ]
  genes <- names(sort(table(gene_list$Name), decreasing = T))
  mat <- matrix(data = 0, nrow = length(genes), ncol = n)
  rownames(mat) <- names(sort(table(gene_list$Name), decreasing = T))
  colnames(mat) <- names
  mat
  for (i in 1:nrow(mat)) {
    print(i)
    for (j in 1:ncol(mat)) {
      mat[i, j] <- any((gene_list$Name == rownames(mat)[i]) & (gene_list$Set ==
                                                                 colnames(mat)[j]))
    }
  }
  
  mat <- data.frame(t(mat))
  write.table(mat,paste("/storage/singlecell/zz4/fetal_snakemake/figures/figure5/",
                        region, "_upset_DEGs.csv",sep=""))
  svg(paste("/storage/singlecell/zz4/fetal_snakemake/figures/figure5/",
            region, "_upset_DEGs.svg",sep=""), width = 15, height = 10, bg = "transparent")
  p <- upset(mat, nintersects = NA, nsets = 15, decreasing = T, order.by = "freq",
             text.scale = 3, keep.order = T, sets = names(sort(table(gene_list$Name),
                                                               decreasing = T)[15:1]), mb.ratio = c(0.4, 0.6), point.size = 5)
  print(p)
  dev.off()
}
mat1 <- read.csv(paste("/storage/singlecell/zz4/fetal_snakemake/figures/figure5/",
                       "Macula", "_upset_DEGs.csv",sep=""),sep = " ")
mat2 <- read.csv(paste("/storage/singlecell/zz4/fetal_snakemake/figures/figure5/",
                       "Peripheral", "_upset_DEGs.csv",sep=""),sep = " ")

data1 <- data.frame(rowSums(mat1))
data1$cell_type <- rownames(data1)
data1$region = "Macula"
data2 <- data.frame(rowSums(mat2))
data2$cell_type <- rownames(data2)
data2$region = "Peripheral"
data <- data.frame(values = c(data1$rowSums.mat1.,data2$rowSums.mat2.),celltype = c(data1$cell_type,data2$cell_type),
                   region = c(data1$region,data2$region))

data$cell_type <- rep(c("PRPC", "PRPC", "PRPC", "NRPC", "NRPC",
                    "NRPC",
                    "MG", "MG",
                    "MG", "AC", "AC", "AC", "BC", "BC",
                    "BC", "Cone", "Cone", "Cone", "Rod", "Rod",
                    "Rod", "HC", "HC", "HC", "RGC", "RGC",
                    "RGC"),2)
data$Time <- rep(c("All","Early", "Late"),9)
data$name <- rownames(data)
data$cell_type <- factor(data$cell_type,levels=c("HC", "MG", "RGC", "AC", "NRPC", "BC", "PRPC", "Cone", "Rod"))
data$celltype <- factor(data$celltype,levels = c(
   "HC_Early", "HC_Late","HC_All",
   "MG_Early","MG_Late","MG_All",
   "RGC_Early","RGC_Late","RGC_All",
   "AC_Early", "AC_Late","AC_All",
   "NRPC_Early","NRPC_Late","NRPC_All",
   "BC_Early","BC_Late","BC_All",
   "PRPC_Early", "PRPC_Late","PRPC_All",
   "Cone_Early", "Cone_Late","Cone_All",
   "Rod_Early","Rod_Late","Rod_All"))
p<-ggplot(data,                         # Draw barplot with grouping & stacking
          aes(x = Time,
              y = values,
              fill = region)) + 
  geom_bar(stat = "identity",
           position = "stack") +
  facet_grid(~ cell_type)+theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ xlab("")+ylab("Number of DEGs")+ guides(fill=guide_legend(title="Region"))
p
