input_path="/storage/chentemp/zz4/adult_dev_compare/data/Retina_fetal/"
output_file="/storage/chentemp/zz4/adult_dev_compare/results/merged_rna/merged_rna.rds"
# Follow vignette at
# https://satijalab.org/seurat/articles/integration_introduction.html
# And https://github.com/satijalab/seurat/issues/1720

# Import packages needed for the following analysis
suppressMessages(library(Seurat))
set.seed(0)
samples <- list.dirs(input_path, full.names = F, recursive = F)

for (sam in samples) {
    # This for loop will read count matrix into Seurat object and
    # Assign Seurat object to each one of sample names.
    cat("\n")
    cat("Processing sample", sam, "now!")
    cat("\n")
    counts <- Read10X_h5(paste0(input_path, sam, "/outs/filtered_feature_bc_matrix.h5",
        sep = ""))
    # Create a Seurat object containing the RNA data
    if (class(counts)=="list") {
    temp <- CreateSeuratObject(counts = counts$`Gene Expression`,
    assay = "RNA")
}
if (class(counts)=="dgCMatrix") {
    temp <- CreateSeuratObject(counts = counts,
    assay = "RNA")
}
    DefaultAssay(temp) <- "RNA"
    temp <- RenameCells(object = temp, add.cell.id = sam)
    assign(paste("A",sam,sep=""), temp)
    cat("Processing sample", sam, "finished!")
    cat("\n")
    cat("\n")
}

seurat_object <- merge(AMulti_Fetal_11w2d_FR, y = c(AMulti_Fetal_11w2d_FR_2, AMulti_Fetal_11w2d_NR, AMulti_Fetal_13W_FR, 
AMulti_Fetal_13W_NR, AMulti_Fetal_14w5d_FR, AMulti_Fetal_14w5d_NR, 
AMulti_Fetal_19W4d_FR, AMulti_Fetal_19W4d_NR, AMulti_Fetal_20W2d_FR, 
AMulti_Fetal_20W2d_NR, AMulti_Fetal_23w1d_FR, AMulti_Fetal_23w1d_NR, 
AMultiome_10w_FR, AMultiome_10w_NR, AMultiome_12w3d_FR, AMultiome_12w3d_NR, 
AMultiome_14w2d_FR, AMultiome_14w2d_NR, AMultiome_16w4d_FR, 
AMultiome_16w4d_NR, AMultiome_20w1d_FR, AMultiome_20w1d_NR, 
AMultiome_23w4d_FR, AMultiome_23w4d_NR, Asn_multiome_d59, 
Asn_multiome_d76c, Asn_multiome_d76p))

saveRDS(seurat_object, output_file)