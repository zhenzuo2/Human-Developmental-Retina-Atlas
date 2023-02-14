subsetbam <- function(bam, cell_id, name, output) {
  sink(paste("/storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_subsetbam/",
             name, ".sh", sep = ""))
  cat(paste("echo \"Processing ", name, "\"\n", sep = ""))
  cat(paste("samtools view -H ", bam, " > ", output, name, "_SAM_header",
            sep = ""))
  cat("\n")
  cat(paste("samtools view ", bam, " | LC_ALL=C grep -F -f ", cell_id,
            " > ", output, name, sep = ""))
  cat("\n")
  cat(paste("cat ", output, name, "_SAM_header ", output, name, " > ",
            output, name, ".sam", sep = ""))
  cat("\n")
  cat(paste("samtools view -b ", output, name, ".sam > ", output, name,
            ".bam", sep = ""))
  cat("\n")
  sink()
}

for (f in list.files("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/earlyPRPC",
                     full.names = T)) {
  print(f)
  name = print(tools::file_path_sans_ext(basename(f)))
  bam = paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
              name, "/outs/atac_possorted_bam.bam", sep = "")
  cell_id = f
  output = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/earlyPRPC/"
  subsetbam(bam, cell_id, name, output)
}

for (f in list.files("/storage/singlecell/zz4/fetal_bash/results/cell_annotation_results/latePRPC",
                     full.names = T)) {
  print(f)
  name = print(tools::file_path_sans_ext(basename(f)))
  bam = paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
              name, "/outs/atac_possorted_bam.bam", sep = "")
  cell_id = f
  output = "/storage/singlecell/zz4/fetal_bash/results/chrombpnet_subsetbam/latePRPC/"
  subsetbam(bam, cell_id, name, output)
}

## cd /storage/singlecell/zz4/fetal_bash/scripts/chrombpnet_subsetbam
## for f in *.sh; do slurmtaco.sh --g01 -m 30G -t 1 -- sh $f; done