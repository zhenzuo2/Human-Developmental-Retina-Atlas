from PIL import Image
import scanpy as sc
import matplotlib.pyplot as plt

def merge_tiff_files(input_files, output_file):
    images = [Image.open(file) for file in input_files]

    # Calculate the number of rows needed
    num_rows = (len(images) + 4) // 5

    # Calculate the size of the output image
    max_width = max(img.width for img in images)
    max_height = max(img.height for img in images)
    output_width = max_width * 5
    output_height = max_height * num_rows

    # Create a new blank image
    output_image = Image.new("RGB", (output_width, output_height), (255, 255, 255))

    # Paste the images into the output image
    for i, img in enumerate(images):
        row = i // 5
        col = i % 5
        x_offset = col * max_width
        y_offset = row * max_height
        output_image.paste(img, (x_offset, y_offset))

    # Save the merged image
    output_image.save(output_file)

early_PRPC = ["SFRP2", "DLX1", "DLX2", "ONECUT2", "ATOH7"]
late_PRPC = ["NFIA", "ASCL1", "OTX2", "SOX4"]
early_NRPC = ["DLX1", "DLX2", "ONECUT2", "ATOH7"]
late_NRPC = ["ASCL1", "OTX2", "SOX4", "OLIG2", "NEUROG2", "BTG2"]
MG = ["SLC1A3", "SLN", "RLBP1", "SOX2", "NFIA", "CRYM", "CLU", "LINC00461"]
adata = sc.read(
    "/storage/chentemp/zz4/adult_dev_compare/results/adata_object_UMAP/ALL.h5ad"
)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

for gene in early_PRPC + late_PRPC + early_NRPC + late_NRPC + MG:
    sc.pl.umap(adata, color=gene, size=10, title=gene,frameon = False)
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.title(gene, fontsize=25)
    plt.savefig(
        "/storage/chentemp/zz4/adult_dev_compare/temp/" + gene + ".tiff",
        bbox_inches="tight",
        transparent=True,
        dpi=300,
    )
# early_PRPC
file_names = [
    "/storage/chentemp/zz4/adult_dev_compare/temp/" + x + ".tiff" for x in early_PRPC
]
output_file_name = (
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/EarlyPRPCMarker.tiff"
)
merge_tiff_files(file_names, output_file_name)

# late_PRPC
file_names = [
    "/storage/chentemp/zz4/adult_dev_compare/temp/" + x + ".tiff" for x in late_PRPC
]
output_file_name = (
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/LatePRPCMarker.tiff"
)
merge_tiff_files(file_names, output_file_name)

# early_NRPC
file_names = [
    "/storage/chentemp/zz4/adult_dev_compare/temp/" + x + ".tiff" for x in early_NRPC
]
output_file_name = (
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/EarlyNRPCMarker.tiff"
)
merge_tiff_files(file_names, output_file_name)

# late_NRPC
file_names = [
    "/storage/chentemp/zz4/adult_dev_compare/temp/" + x + ".tiff" for x in late_NRPC
]
output_file_name = (
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/late_NRPCMarker.tiff"
)
merge_tiff_files(file_names, output_file_name)

# MG
file_names = ["/storage/chentemp/zz4/adult_dev_compare/temp/" + x + ".tiff" for x in MG]
output_file_name = (
    "/storage/chentemp/zz4/adult_dev_compare/Supplementary/Figure1/MGMarker.tiff"
)
merge_tiff_files(file_names, output_file_name)
