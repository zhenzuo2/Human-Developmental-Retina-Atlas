mapping <- read.csv("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal_sample_meta.csv")
mapping <- rbind(mapping, c("Multi_Fetal_11w2d_FRre", "11w2d", "Macula",
                            79))
mapping <- rbind(mapping, c("Multi_organoid_Fetal_11w2d_NR", "11w2d", "Peripheral",
                            79))
mapping <- rbind(mapping, c("sample_sheet_Multi_Fetal_14w5d_NR", "14w5d",
                            "Peripheral", 103))
mapping <- rbind(mapping, c("Multi_12w3d_FR", "12w2d", "Macula", 87))
mapping <- rbind(mapping, c("Multi_12w3d_NR", "12w3d", "Peripheral", 87))
mapping <- rbind(mapping, c("10xMulti_23w4d", "23w4d", "Peripheral", 165))
samples <- list.dirs("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                     full.names = F, recursive = F)
samples = samples[4:28]
df = data.frame()
for (sample in samples) {
  print(sample)
  if (file.exists(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                        sample, "/outs/metrics_summary.csv", sep = ""))) {
    df <- plyr::rbind.fill(df, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/metrics_summary.csv", sep = "")))
  } else {
    df <- plyr::rbind.fill(df, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/summary.csv", sep = "")))
  }
}
df <- df[, 1:46]
df$Time <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Time)
df$Region <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Region)
df$Days <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Days)

df <- df[, c("Sample.ID", "Genome", "Pipeline.version", "Time", "Days",
             "Region", "Estimated.number.of.cells", "Feature.linkages.detected",
             "Linked.genes", "Linked.peaks", "ATAC.Confidently.mapped.read.pairs",
             "ATAC.Fraction.of.genome.in.peaks", "ATAC.Fraction.of.high.quality.fragments.in.cells",
             "ATAC.Fraction.of.high.quality.fragments.overlapping.TSS", "ATAC.Fraction.of.high.quality.fragments.overlapping.peaks",
             "ATAC.Fraction.of.transposition.events.in.peaks.in.cells", "ATAC.Mean.raw.read.pairs.per.cell",
             "ATAC.Median.high.quality.fragments.per.cell", "ATAC.Non.nuclear.read.pairs",
             "ATAC.Number.of.peaks", "ATAC.Percent.duplicates", "ATAC.Q30.bases.in.barcode",
             "ATAC.Q30.bases.in.read.1", "ATAC.Q30.bases.in.read.2", "ATAC.Q30.bases.in.sample.index.i1",
             "ATAC.Sequenced.read.pairs", "ATAC.TSS.enrichment.score", "ATAC.Unmapped.read.pairs",
             "ATAC.Valid.barcodes", "GEX.Fraction.of.transcriptomic.reads.in.cells",
             "GEX.Mean.raw.reads.per.cell", "GEX.Median.UMI.counts.per.cell", "GEX.Median.genes.per.cell",
             "GEX.Percent.duplicates", "GEX.Q30.bases.in.UMI", "GEX.Q30.bases.in.barcode",
             "GEX.Q30.bases.in.read.2", "GEX.Reads.mapped.antisense.to.gene", "GEX.Reads.mapped.confidently.to.exonic.regions",
             "GEX.Reads.mapped.confidently.to.genome", "GEX.Reads.mapped.confidently.to.intergenic.regions",
             "GEX.Reads.mapped.confidently.to.intronic.regions", "GEX.Reads.mapped.confidently.to.transcriptome",
             "GEX.Reads.mapped.to.genome", "GEX.Reads.with.TSO", "GEX.Sequenced.read.pairs",
             "GEX.Total.genes.detected", "GEX.Valid.UMIs", "GEX.Valid.barcodes")]

write.table(df,"/storage/singlecell/zz4/fetal_bash/results/meta/first_round_meta.csv")

samples <- list.dirs("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                     full.names = F, recursive = F)
samples = samples[1:3]
df = data.frame()
for (sample in samples) {
  print(sample)
  if (file.exists(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                        sample, "/outs/metrics_summary.csv", sep = ""))) {
    df <- plyr::rbind.fill(df, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/metrics_summary.csv", sep = "")))
  } else {
    df <- plyr::rbind.fill(df, read.csv(paste("/storage/singlecell/zz4/fetal_bash/data/Retina_fetal/",
                                              sample, "/outs/summary.csv", sep = "")))
  }
}
df$Sample.ID = samples
df$Time <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Time)
df$Region <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Region)
df$Days <- plyr::mapvalues(df$Sample.ID, from = mapping$Samples, to = mapping$Days)
df <- df[, c("Sample.ID", "Time", "Region", "Days", "Estimated.Number.of.Cells",
             "Mean.Reads.per.Cell", "Median.Genes.per.Cell", "Number.of.Reads",
             "Valid.Barcodes", "Sequencing.Saturation", "Q30.Bases.in.Barcode",
             "Q30.Bases.in.RNA.Read", "Q30.Bases.in.UMI", "Reads.Mapped.to.Genome",
             "Reads.Mapped.Confidently.to.Genome", "Reads.Mapped.Confidently.to.Intergenic.Regions",
             "Reads.Mapped.Confidently.to.Intronic.Regions", "Reads.Mapped.Confidently.to.Exonic.Regions",
             "Reads.Mapped.Confidently.to.Transcriptome", "Reads.Mapped.Antisense.to.Gene",
             "Fraction.Reads.in.Cells", "Total.Genes.Detected", "Median.UMI.Counts.per.Cell")]

write.table(df,"/storage/singlecell/zz4/fetal_bash/results/meta/second_round.csv")