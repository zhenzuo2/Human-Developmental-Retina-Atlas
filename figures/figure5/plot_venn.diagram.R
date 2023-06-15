library(ggplot2)
library(RColorBrewer)
library(VennDiagram)

output_dir = c("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_",
  "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_",
  "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_")

cm_files = c("All_compare_mod_region.csv", "Early_compare_mod_region.csv",
             "Late_compare_mod_region.csv")

coef_files = c("All_fit_coefs_region_models.csv", "Early_fit_coefs_region_models.csv",
               "Late_fit_coefs_region_models.csv")

mode = c("All", "Early", "Late")

for (out in output_dir) {
  for (i in 1:3) {
    coef <- read.csv(paste(out, coef_files[i], sep = ""))
    coef <- coef[(coef$status == "OK") & (coef$term == "RegionPeripheral") &
                   (coef$q_value < 0.01) & (abs(coef$normalized_effect) >=0.5),
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

files <- c("/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/PRPC_Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/NRPC_Late_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_All_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_Early_Region_DE_filtered_gene_list.csv",
           "/storage/singlecell/zz4/fetal_snakemake/results/monocle3_DE_analysis/MG_monocle3_DE_analysis/RPC_Late_Region_DE_filtered_gene_list.csv")

names <- c("PRPC_All", "PRPC_Early",
           "PRPC_Late", "NRPC_All", "NRPC_Early", "NRPC_Late", "RPC_All", "RPC_Early", "RPC_Late")

n <- length(files)
gene_list <- data.frame(Name = character(), Set = character())
for (i in 1:n) {
  temp <- read.csv(files[i])
  colnames(temp) <- "Name"
  temp$Set <- names[i]
  gene_list <- rbind(gene_list, temp)
}
gene_list <- gene_list[!grepl("[^A-Za-z0-9]", gene_list$Name),]
genes <- names(sort(table(gene_list$Name),decreasing = T))

# Generate 3 sets of 200 words
set1 <- gene_list[gene_list$Set%in%c("PRPC_Early"),"Name"]
set2 <- gene_list[gene_list$Set%in%c("PRPC_Late"),"Name"]
set3 <- gene_list[gene_list$Set%in%c("NRPC_Early"),"Name"]
set4 <- gene_list[gene_list$Set%in%c("NRPC_Late"),"Name"]


# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
# Load library
library(VennDiagram)

# Chart
venn.diagram(
  x = list(set1, set2,set3,set4),
  category.names = c("Early PRPC" , "Late PRPC","Early NRPC" , "Late NRPC"),
  filename = '/storage/singlecell/zz4/fetal_snakemake/figures/figure5/venn_diagramm.png',
  fill = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  output=TRUE,cex = 1.5,cat.cex = 1.5,margin=0.1
)

