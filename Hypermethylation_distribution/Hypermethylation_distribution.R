library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPseeker)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# this file has 3 columns: sample, bigwig, NAPY
# sample column will be the name sued during plotting
# bigwig column is the path to the bigwig file
# NAPY: is an annotation column with the NAPY status
cell_line_bws = read.table("Cell_line_bigwigs_withNAPY.tsv", header = T, sep = "\t")
cell_line_bws = cell_line_bws[cell_line_bws$NAPY != "NSCLC",]

# Importing bigwig files
bws = list()
for(i in 1:nrow(cell_line_bws)) {
  print(paste0("Importing ", i, " of ", nrow(cell_line_bws)))
  bws[[i]] = rtracklayer::import(cell_line_bws$bigwig[i])
}

# Combining bigwig files into a data frame
probe.df = data.frame(matrix(ncol =  nrow(cell_line_bws), nrow = length(bws[[1]])))
refprobe = bws[[1]]
names(refprobe) = seq(1:length(refprobe))
colnames(probe.df) =  cell_line_bws$sample
rownames(probe.df) = names(refprobe)
for(i in 1:nrow(cell_line_bws)) {
  bw = bws[[i]]
  overs = findOverlaps(query = bws[[i]], subject = refprobe)
  probe.df[subjectHits(overs),i] = bws[[i]][queryHits(overs)]$score
}
probe.df = probe.df[complete.cases(probe.df),]
refprobe = refprobe[rownames(probe.df),]

# Annotation of probes
refprobe$annotation = "Intergenic"

# getting gene annotating probes
overs = findOverlaps(query = refprobe, subject = genes(TxDb.Hsapiens.UCSC.hg19.knownGene))
refprobe[queryHits(overs)]$annotation = "Gene body"

# getting promoter annotating probes
overs = findOverlaps(query = refprobe, subject = promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = promoter_extend_bases, downstream = promoter_extend_bases))
refprobe[queryHits(overs)]$annotation = "Promoter"

# getting CpG annotating probes
overs = findOverlaps(query = refprobe, subject = cpgs)
refprobe[queryHits(overs)]$annotation = "CpG"


plot.df = reshape2::melt(as.matrix(probe.df))
plot.df$Var1 = as.character(plot.df$Var1)
plot.df$Var2 = as.character(plot.df$Var2)
plot.df$annotation = refprobe[as.character(plot.df$Var1)]$annotation
plot.df$methylation_status = ifelse(plot.df$value > 0.5, "high", "low")

cell_table = table(plot.df[,c("Var2","methylation_status", "annotation")])

cell_lines = as.character(unique(plot.df$Var2))
annotation_types = as.character(unique(plot.df$annotation))

anndf = data.frame()

for(i in 1:length(cell_lines)) {
  print(paste0("Analyzing ", i, " of ", nrow(cell_line_bws)))
  for(j in 1:length(annotation_types)) {
    tmp = plot.df[plot.df$Var2 == cell_lines[i] & plot.df$annotation == annotation_types[j],]
    
    if(nrow(anndf) == 0) {
      anndf = data.frame(cell = cell_lines[i], annotation = annotation_types[j], 
                         hyper_fraction = 100 * nrow(tmp[tmp$methylation_status == "high",])/ nrow(tmp))
    } else {
      anndf = rbind(anndf, data.frame(cell = cell_lines[i], annotation = annotation_types[j], 
                                      hyper_fraction = 100 * nrow(tmp[tmp$methylation_status == "high",])/ nrow(tmp)))
    }
  }
}

for(i in 1:length(cell_lines)) {
  tmp = plot.df[plot.df$Var2 == cell_lines[i],]
  anndf = rbind(anndf, data.frame(cell = cell_lines[i], annotation = "All", 
                                  hyper_fraction = 100 * nrow(tmp[tmp$methylation_status == "high",])/ nrow(tmp)))
}

write.table(anndf, "Hypermethylation_fractions.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
