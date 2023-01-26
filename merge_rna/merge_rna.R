args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_file <- args[2]
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
    assign(sam, temp)
    cat("Processing sample", sam, "finished!")
    cat("\n")
    cat("\n")
}

seurat_object <- merge(get("17W1D_Fovea_retina"), y = c(get("17w1d_I_Ret"), get("17W1D_Nasal_retina"), 
get("17w1d_S_Ret"), get("17W1D_Temporal_retina"), multi_19w3d_F_ret, 
multi_19W3D_I_ret, multi_19W3D_N_RET, multi_19w3d_S_ret, 
multi_19W3d_T_ret, Multi_Fetal_11w2d_FR, Multi_Fetal_11w2d_FR_2, 
Multi_Fetal_11w2d_NR, Multi_Fetal_13W_FR, Multi_Fetal_13W_NR, 
Multi_Fetal_14w5d_FR, Multi_Fetal_14w5d_NR, Multi_Fetal_19W4d_FR, 
Multi_Fetal_19W4d_NR, Multi_Fetal_20W2d_FR, Multi_Fetal_20W2d_NR, 
Multi_Fetal_23w1d_FR, Multi_Fetal_23w1d_NR, Multiome_10w_FR, 
Multiome_10w_NR, Multiome_12w3d_FR, Multiome_12w3d_NR, 
Multiome_14w2d_FR, Multiome_14w2d_NR, Multiome_16w4d_FR, 
Multiome_16w4d_NR, Multiome_20w1d_FR, Multiome_20w1d_NR, 
Multiome_23w4d_FR, Multiome_23w4d_NR))

saveRDS(seurat_object, output_file)